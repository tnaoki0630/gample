#import "EMField.h"
#import <Metal/Metal.h>

@interface EMField () {
    id<MTLDevice> _device;
    float* _rho;                 // 電荷密度
    float* _phi;                 // ポテンシャル
    id<MTLBuffer> _ExBuffer;     // 電場 x成分
    id<MTLBuffer> _EyBuffer;     // 電場 y成分
    id<MTLBuffer> _EzBuffer;     // 電場 z成分
    id<MTLBuffer> _BxBuffer;     // 磁場 x成分
    id<MTLBuffer> _ByBuffer;     // 磁場 y成分
    id<MTLBuffer> _BzBuffer;     // 磁場 z成分
    
    int _ngx;                    // x方向グリッド数
    int _ngy;                    // y方向グリッド数
    int _ngb;                    // バッファ領域グリッド数
    float _dx;                  // x方向グリッド幅
    float _dy;                  // y方向グリッド幅
    
    int _weightOrder;            // weighting のオーダー
}
@end

@implementation EMField {
    // プライベートインスタンス変数
    int _ngbrho[2];
    int _ngbphi[2];
    int _ngbEx[2];
    int _ngbEy[2];
    int _ngbEz[2];
    int _ngbBx[2];
    int _ngbBy[2];
    int _ngbBz[2];
}

- (instancetype)initWithDevice:(id<MTLDevice>)device withParam:(Init*)initParam {
    self = [super init];
    if (self) {
        _device = device;
        
        // フィールドパラメータの取得
        struct ParamForField fieldParam = [initParam getParamForField];
        _ngx = fieldParam.ngx;
        _ngy = fieldParam.ngy;
        _dx = fieldParam.dx;
        _dy = fieldParam.dy;
        _weightOrder = fieldParam.weightOrder;

        // 小行列の配列サイズ（A[1][1] が解析領域で、それ以外は領域境界）
        if (_weightOrder == 5){
            // general
            _ngb = 2;
            // rho,phi,Ez
            _ngbrho[0] = 2;
            _ngbphi[0] = 2;
            _ngbEz[0] = 2;
            _ngbrho[1] = 2;
            _ngbphi[1] = 2;
            _ngbEz[1] = 2;
            // Ex, By
            _ngbEx[0] = 3;
            _ngbBy[0] = 3;
            _ngbEx[1] = 2;
            _ngbBy[1] = 2;
            // Ey, Bx
            _ngbEy[0] = 2;
            _ngbBx[0] = 2;
            _ngbEy[1] = 3;
            _ngbBx[1] = 3;
            // Bz
            _ngbBz[0] = 3;
            _ngbBz[0] = 3;
        } else {
            _ngb = 0;
        }

        // バッファサイズの計算
        NSUInteger gridSize = (_ngx + 2*_ngb) * (_ngy + 2*_ngb) ;
        NSUInteger bufferSize = sizeof(float) * gridSize;
        
        // Metal バッファの初期化
        _ExBuffer = [device newBufferWithLength:bufferSize options:MTLResourceStorageModeShared];
        _EyBuffer = [device newBufferWithLength:bufferSize options:MTLResourceStorageModeShared];
        _EzBuffer = [device newBufferWithLength:bufferSize options:MTLResourceStorageModeShared];
        _BxBuffer = [device newBufferWithLength:bufferSize options:MTLResourceStorageModeShared];
        _ByBuffer = [device newBufferWithLength:bufferSize options:MTLResourceStorageModeShared];
        _BzBuffer = [device newBufferWithLength:bufferSize options:MTLResourceStorageModeShared];
        
        // バッファの初期化
        float* Ex = (float*)_ExBuffer.contents;
        float* Ey = (float*)_EyBuffer.contents;
        float* Ez = (float*)_EzBuffer.contents;
        float* Bx = (float*)_BxBuffer.contents;
        float* By = (float*)_ByBuffer.contents;
        float* Bz = (float*)_BzBuffer.contents;
        
        // malloc
        _rho = (float *)malloc(bufferSize);
        _phi = (float *)malloc(bufferSize);
        
        // initialize(Uniform)
        for (NSUInteger i = 0; i < gridSize; i++) {
            _rho[i] = 0.0f;
            _phi[i] = 0.0f;
            Ex[i] = fieldParam.ampE[0];
            Ey[i] = fieldParam.ampE[1];
            Ez[i] = fieldParam.ampE[2];
            Bx[i] = fieldParam.ampB[0];
            By[i] = fieldParam.ampB[1];
            Bz[i] = fieldParam.ampB[2];
        }
    }
    return self;
}

- (void)solvePoisson{

}

- (void)resetChargeDensity{
    // 電荷密度の初期化
    int ng = (_ngx + 2*_ngb) * (_ngy + 2*_ngb);
    for (int i = 0; i < ng; i++){
        _rho[i] = 0.0f;
    }
}

- (void)outputField:(int)cycle{
    // サイクル番号を用いてファイル名を作成 (例: field_0001.bin)
    NSString *fileName = [NSString stringWithFormat:@"bin/field_%08d.bin", cycle];
    const char *filePath = [fileName UTF8String];
    
    FILE *fp = fopen(filePath, "wb");
    if (!fp) {
        NSLog(@"Error: Unable to open file %s for writing", filePath);
        return;
    }
    
    // ヘッダ情報としてメッシュ情報を出力
    fwrite(&_ngx, sizeof(int), 1, fp);
    fwrite(&_ngy, sizeof(int), 1, fp);
    fwrite(&_ngb, sizeof(int), 1, fp);
    fwrite(&_dx, sizeof(float), 1, fp);
    fwrite(&_dy, sizeof(float), 1, fp);
    
    // 全グリッドサイズ (境界含む) の計算
    int ng = (_ngx + 2 * _ngb) * (_ngy + 2 * _ngb);
    
    // 順にフィールドデータを書き出す: _rho, _phi, Ex, Ey, Bz
    fwrite(_rho, sizeof(float), ng, fp);
    fwrite(_phi, sizeof(float), ng, fp);
    
    float *Ex = (float *)_ExBuffer.contents;
    float *Ey = (float *)_EyBuffer.contents;
    float *Bz = (float *)_BzBuffer.contents;
    
    fwrite(Ex, sizeof(float), ng, fp);
    fwrite(Ey, sizeof(float), ng, fp);
    fwrite(Bz, sizeof(float), ng, fp);
    
    fclose(fp);
    NSLog(@"Field data successfully written to %s", filePath);
}

// 電荷密度へのアクセサ
- (float*)rho { return _rho; }

// 電場バッファへのアクセサ
- (id<MTLBuffer>)ExBuffer { return _ExBuffer; }
- (id<MTLBuffer>)EyBuffer { return _EyBuffer; }
- (id<MTLBuffer>)EzBuffer { return _EzBuffer; }

// 磁場バッファへのアクセサ
- (id<MTLBuffer>)BxBuffer { return _BxBuffer; }
- (id<MTLBuffer>)ByBuffer { return _ByBuffer; }
- (id<MTLBuffer>)BzBuffer { return _BzBuffer; }

// グリッド情報へのアクセサ
- (int)ngx { return _ngx; }
- (int)ngy { return _ngy; }
- (int)ngb { return _ngb; }
- (double)dx { return _dx; }
- (double)dy { return _dy; }

@end

#import "EMField.h"