#import "PIC.h"

static NSString *const kMetalShaderSource = @R"(
#include <metal_stdlib>
using namespace metal;

struct Particle {
    float2 position;
    float2 velocity;
    float charge;
    float mass;
};

struct Field {
    float2 value;
};

// 粒子更新カーネル
kernel void updateParticles(device Particle* particles [[buffer(0)]],
                          device const Field* fields [[buffer(1)]],
                          const device uint& particleCount [[buffer(2)]],
                          const device float& deltaTime [[buffer(3)]],
                          const device float2& gridSize [[buffer(4)]],
                          uint id [[thread_position_in_grid]]) {
    if (id >= particleCount) return;
    
    // 参照ではなくポインタとして扱う
    device Particle* p = &particles[id];
    
    // 最近接格子点からの電場を計算
    uint2 gridPos = uint2(p->position);
    Field localField = fields[gridPos.y * uint(gridSize.x) + gridPos.x];
    
    // 運動方程式の解析
    float2 acceleration = (p->charge / p->mass) * localField.value;
    p->velocity += acceleration * deltaTime;
    p->position += p->velocity * deltaTime;
    
    // 周期的境界条件
    p->position = fmod(p->position + gridSize, gridSize);
}

// 電場計算カーネル
kernel void computeFields(device Field* fields [[buffer(0)]],
                        device const Particle* particles [[buffer(1)]],
                        const device uint& particleCount [[buffer(2)]],
                        const device float2& gridSize [[buffer(3)]],
                        uint2 pos [[thread_position_in_grid]]) {
    if (pos.x >= uint(gridSize.x) || pos.y >= uint(gridSize.y)) return;
    
    uint index = pos.y * uint(gridSize.x) + pos.x;
    float2 fieldValue = float2(0, 0);
    
    // 各粒子からのクーロン力を計算
    for (uint i = 0; i < particleCount; i++) {
        // 参照ではなくポインタとして直接アクセス
        const device Particle* p = &particles[i];
        float2 r = p->position - float2(pos);
        float r2 = dot(r, r);
        if (r2 > 0.0001f) {  // 特異点を避ける
            fieldValue += p->charge * normalize(r) / r2;
        }
    }
    
    fields[index].value = fieldValue;
}
)";

@implementation PIC

- (instancetype)initWithDevice:(id<MTLDevice>)device
                particleCount:(NSUInteger)particleCount
                    gridSize:(NSUInteger)gridSize {
    self = [super init];
    if (self) {
        _device = device;
        _commandQueue = [device newCommandQueue];
        
        // コンピュートパイプラインの設定
        NSError *error = nil;
        id<MTLLibrary> library = [device newLibraryWithSource:kMetalShaderSource
                                                     options:nil
                                                       error:&error];
        if (!library) {
            NSLog(@"Failed to create Metal library: %@", error);
            return nil;
        }
        
        id<MTLFunction> updateParticlesFunction = [library newFunctionWithName:@"updateParticles"];
        id<MTLFunction> computeFieldsFunction = [library newFunctionWithName:@"computeFields"];
        
        _updateParticlesPipeline = [device newComputePipelineStateWithFunction:updateParticlesFunction
                                                                        error:&error];
        _computeFieldsPipeline = [device newComputePipelineStateWithFunction:computeFieldsFunction
                                                                      error:&error];
        
        // バッファの初期化
        _particleBuffer = [device newBufferWithLength:sizeof(Particle) * particleCount
                                            options:MTLResourceStorageModeShared];
        _fieldBuffer = [device newBufferWithLength:sizeof(Field) * gridSize * gridSize
                                         options:MTLResourceStorageModeShared];
        
        // 初期粒子分布の設定
        [self initializeParticles:particleCount gridSize:gridSize];
    }
    return self;
}

- (void)initializeParticles:(NSUInteger)particleCount gridSize:(NSUInteger)gridSize {
    Particle *particles = (Particle *)self.particleBuffer.contents;
    
    for (NSUInteger i = 0; i < particleCount; i++) {
        particles[i].position = (vector_float2){
            (float)(arc4random_uniform(gridSize)),
            (float)(arc4random_uniform(gridSize))
        };
        particles[i].velocity = (vector_float2){0.0f, 0.0f};
        particles[i].charge = (i % 2) * 2.0f - 1.0f;  // 正負交互に配置
        particles[i].mass = 1.0f;
    }
}

- (void)update:(float)deltaTime {
    id<MTLCommandBuffer> commandBuffer = [self.commandQueue commandBuffer];
    
    // 電場の計算
    {
        id<MTLComputeCommandEncoder> computeEncoder = [commandBuffer computeCommandEncoder];
        [computeEncoder setComputePipelineState:self.computeFieldsPipeline];
        [computeEncoder setBuffer:self.fieldBuffer offset:0 atIndex:0];
        [computeEncoder setBuffer:self.particleBuffer offset:0 atIndex:1];
        // ... その他のパラメータの設定 ...
        
        MTLSize gridSize = MTLSizeMake(32, 32, 1);
        MTLSize threadGroupSize = MTLSizeMake(8, 8, 1);
        [computeEncoder dispatchThreadgroups:gridSize
                    threadsPerThreadgroup:threadGroupSize];
        [computeEncoder endEncoding];
    }
    
    // 粒子の更新
    {
        id<MTLComputeCommandEncoder> computeEncoder = [commandBuffer computeCommandEncoder];
        [computeEncoder setComputePipelineState:self.updateParticlesPipeline];
        [computeEncoder setBuffer:self.particleBuffer offset:0 atIndex:0];
        [computeEncoder setBuffer:self.fieldBuffer offset:0 atIndex:1];
        // ... その他のパラメータの設定 ...
        
        MTLSize gridSize = MTLSizeMake(1024, 1, 1);
        MTLSize threadGroupSize = MTLSizeMake(256, 1, 1);
        [computeEncoder dispatchThreadgroups:gridSize
                    threadsPerThreadgroup:threadGroupSize];
        [computeEncoder endEncoding];
    }
    
    [commandBuffer commit];
}

- (void)writeVTKFile:(NSString *)filename forTimestep:(NSInteger)timestep {
    // ファイル名の生成
    NSString *fullPath = [filename stringByAppendingFormat:@"_%04ld.vtp", (long)timestep];
    const char *cPath = [fullPath UTF8String];
    
    // ファイルを開く
    std::ofstream file(cPath);
    if (!file.is_open()) {
        NSLog(@"Failed to open file: %@", fullPath);
        return;
    }
    
    // パーティクルデータへのアクセス
    Particle *particles = (Particle *)self.particleBuffer.contents;
    NSUInteger particleCount = self.particleBuffer.length / sizeof(Particle);
    
    // VTK XML ヘッダー
    file << "<?xml version=\"1.0\"?>\n";
    file << "<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    file << "  <PolyData>\n";
    file << "    <Piece NumberOfPoints=\"" << particleCount << "\" NumberOfVerts=\"" << particleCount << "\">\n";
    
    // 座標データ
    file << "      <Points>\n";
    file << "        <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (NSUInteger i = 0; i < particleCount; i++) {
        file << "          " << particles[i].position.x << " " 
             << particles[i].position.y << " 0.0\n";
    }
    file << "        </DataArray>\n";
    file << "      </Points>\n";
    
    // 頂点データ
    file << "      <Verts>\n";
    file << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    for (NSUInteger i = 0; i < particleCount; i++) {
        file << "          " << i << "\n";
    }
    file << "        </DataArray>\n";
    file << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    for (NSUInteger i = 0; i < particleCount; i++) {
        file << "          " << i + 1 << "\n";
    }
    file << "        </DataArray>\n";
    file << "      </Verts>\n";
    
    // 物理量データ
    file << "      <PointData>\n";
    // 速度
    file << "        <DataArray type=\"Float32\" Name=\"velocity\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (NSUInteger i = 0; i < particleCount; i++) {
        file << "          " << particles[i].velocity.x << " " 
             << particles[i].velocity.y << " 0.0\n";
    }
    file << "        </DataArray>\n";
    // 電荷
    file << "        <DataArray type=\"Float32\" Name=\"charge\" format=\"ascii\">\n";
    for (NSUInteger i = 0; i < particleCount; i++) {
        file << "          " << particles[i].charge << "\n";
    }
    file << "        </DataArray>\n";
    file << "      </PointData>\n";
    
    // フッター
    file << "    </Piece>\n";
    file << "  </PolyData>\n";
    file << "</VTKFile>\n";
    
    file.close();
}

- (void)writeFieldVTKFile:(NSString *)filename forTimestep:(NSInteger)timestep {
    NSString *fullPath = [filename stringByAppendingFormat:@"_field_%04ld.vti", (long)timestep];
    const char *cPath = [fullPath UTF8String];
    
    std::ofstream file(cPath);
    if (!file.is_open()) {
        NSLog(@"Failed to open field file: %@", fullPath);
        return;
    }
    
    // フィールドデータへのアクセス
    Field *fields = (Field *)self.fieldBuffer.contents;
    NSUInteger gridSize = (NSUInteger)sqrt(self.fieldBuffer.length / sizeof(Field));
    
    // VTK XML ヘッダー
    file << "<?xml version=\"1.0\"?>\n";
    file << "<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    file << "  <ImageData WholeExtent=\"0 " << gridSize-1 << " 0 " << gridSize-1 << " 0 0\"\n";
    file << "             Origin=\"0 0 0\"\n";
    file << "             Spacing=\"1 1 1\">\n";
    file << "    <Piece Extent=\"0 " << gridSize-1 << " 0 " << gridSize-1 << " 0 0\">\n";
    
    // 電場ベクトルデータ
    file << "      <PointData>\n";
    file << "        <DataArray type=\"Float32\" Name=\"E\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (NSUInteger y = 0; y < gridSize; y++) {
        for (NSUInteger x = 0; x < gridSize; x++) {
            NSUInteger index = y * gridSize + x;
            file << "          " << fields[index].value.x << " "
                 << fields[index].value.y << " 0.0\n";
        }
    }
    file << "        </DataArray>\n";
    
    // 電場強度（スカラー場）の出力
    file << "        <DataArray type=\"Float32\" Name=\"E_magnitude\" format=\"ascii\">\n";
    for (NSUInteger y = 0; y < gridSize; y++) {
        for (NSUInteger x = 0; x < gridSize; x++) {
            NSUInteger index = y * gridSize + x;
            float magnitude = sqrt(fields[index].value.x * fields[index].value.x +
                                 fields[index].value.y * fields[index].value.y);
            file << "          " << magnitude << "\n";
        }
    }
    file << "        </DataArray>\n";
    file << "      </PointData>\n";
    
    // フッター
    file << "    </Piece>\n";
    file << "  </ImageData>\n";
    file << "</VTKFile>\n";
    
    file.close();
}

@end
