#import "FileIO.h"
#import "Constant.h"

void outputPhaseSpace(int cycle, Particle* ptcl, Init* init, XmlLogger& logger){
    struct ParamForField fldParam = init.paramForField;
    struct ParamForTimeIntegration timeParam = init.paramForTimeIntegration;
    
    // obtain objects
    ParticleState *p = (ParticleState*)[[ptcl particleBuffer] contents];
    SimulationParams *prm = (SimulationParams*)[[ptcl paramsBuffer] contents];

    float x,y;
    // 各粒子についてバイナリファイルに出力する 
    for (int idx = 0; idx < 20; idx++) {
        NSString *filePath = [NSString stringWithFormat:@"bin/%@_PhaseSpace_%@_%d.bin", timeParam.ProjectName, ptcl.pName, idx];
        // バイナリ書き出し
        std::ofstream ofs([filePath UTF8String], std::ios::binary | std::ios::app);
        if (!ofs) {
            NSLog(@"Failed to open file: %@", filePath);
            continue;
        }

        // 位置を物理次元に戻す
        x = (p[idx].x - float(fldParam.ngb))*fldParam.dx;
        y = (p[idx].y - float(fldParam.ngb))*fldParam.dy;
        
        // phasespace を出力
        ofs.write(reinterpret_cast<const char*>(&cycle), sizeof(int));
        ofs.write(reinterpret_cast<const char*>(&x), sizeof(float));
        ofs.write(reinterpret_cast<const char*>(&y), sizeof(float));
        ofs.write(reinterpret_cast<const char*>(&p[idx].vx), sizeof(float));
        ofs.write(reinterpret_cast<const char*>(&p[idx].vy), sizeof(float));
        ofs.write(reinterpret_cast<const char*>(&p[idx].vz), sizeof(float));
        
        ofs.close();
        NSLog(@"PhaseSpace successfully written to %@", filePath);
    }

}

void outputField(int cycle, EMField* fld, Init* init, XmlLogger& logger){
    struct ParamForTimeIntegration timeParam = init.paramForTimeIntegration;

    NSString *fileName = [NSString stringWithFormat:@"bin/%@_EMField_%08d.bin", timeParam.ProjectName, cycle];
    const char *filePath = [fileName UTF8String];

    // fld クラスの内容を取得
    float* phi = fld.phi;
    float* rho = (float *)[fld.rhoBuffer contents];
    float* Ex  = (float *)[fld.ExBuffer  contents]; 
    float* Ey  = (float *)[fld.EyBuffer  contents]; 
    float* Bz  = (float *)[fld.BzBuffer  contents]; 

    // グリッド情報取得
    int ngx = fld.ngx;
    int ngy = fld.ngy;
    int ngb = fld.ngb;
    int nx = ngx+2*ngb;
    int ny = ngy+2*ngb;
    float dx = fld.dx;
    float dy = fld.dy;

    // minmax
    float min_rho = 1e20, max_rho = -1e20;
    float min_phi = 1e20, max_phi = -1e20;
    float min_Ex = 1e20, max_Ex = -1e20;
    float min_Ey = 1e20, max_Ey = -1e20;
    for(int i = 0; i < (nx+1)*(ny+1); i++){
        if (rho[i] < min_rho)       { min_rho = rho[i]; }
        else if(rho[i] > max_rho)   { max_rho = rho[i]; }
    }
    for(int i = 0; i < (nx+2)*(ny+1); i++){
        if (Ex[i] < min_Ex)         { min_Ex = Ex[i]; }
        else if(Ex[i] > max_Ex)     { max_Ex = Ex[i]; }
    }
    for(int i = 0; i < (nx+1)*(ny+2); i++){
        if (Ey[i] < min_Ey)         { min_Ey = Ey[i]; }
        else if(Ey[i] > max_Ey)     { max_Ey = Ey[i]; }
    }
    for(int i = 0; i < (nx+3)*(ny+3); i++){
        if (phi[i] < min_phi)       { min_phi = phi[i]; }
        else if(phi[i] > max_phi)   { max_phi = phi[i]; }
    }
    // NSLog(@"MinMax: min_rho = %e, max_rho = %e", min_rho, max_rho);
    std::map<std::string, std::string> data ={
        {"rho_min", fmtSci(min_rho, 6)},
        {"rho_max", fmtSci(max_rho, 6)},
        {"phi_min", fmtSci(min_phi*sVtoV, 6)},
        {"phi_max", fmtSci(max_phi*sVtoV, 6)},
        {"Ex_min", fmtSci(min_Ex*GtoV, 6)},
        {"Ex_max", fmtSci(max_Ex*GtoV, 6)},
        {"Ey_min", fmtSci(min_Ey*GtoV, 6)},
        {"Ey_max", fmtSci(max_Ey*GtoV, 6)},
    };
    logger.logSection("outputField", data);
    
    FILE *fp = fopen(filePath, "wb");
    if (!fp) {
        NSLog(@"Error: Unable to open file %s for writing", filePath);
        return;
    }
    
    // ヘッダ情報としてメッシュ情報を出力
    fwrite(&ngx, sizeof(int), 1, fp);
    fwrite(&ngy, sizeof(int), 1, fp);
    fwrite(&ngb, sizeof(int), 1, fp);
    fwrite(&dx , sizeof(float), 1, fp);
    fwrite(&dy , sizeof(float), 1, fp);
    
    // フィールドデータを書き出す: name,type,array
    writeField(fp, "rho", 0, rho, (nx+1)*(ny+1), 1.0f);
    writeField(fp, "phi", 4, phi, (nx+3)*(ny+3), sVtoV);
    writeField(fp, "Ex", 1, Ex, (nx+2)*(ny+1), GtoV);
    writeField(fp, "Ey", 2, Ey, (nx+1)*(ny+2), GtoV);
    writeField(fp, "Bz", 3, Bz, (nx+2)*(ny+2), GtoT);
    
    fclose(fp);
    NSLog(@"Field data successfully written to %s", filePath);
}

void outputMoments(int cycle, Moment* mom, NSString* pName, Init* init, XmlLogger& logger){
    struct ParamForField fldParam = init.paramForField;
    struct ParamForTimeIntegration timeParam = init.paramForTimeIntegration;

    NSString *fileName = [NSString stringWithFormat:@"bin/%@_Moments_%@_%08d.bin", timeParam.ProjectName, pName, cycle];
    const char *filePath = [fileName UTF8String];

    // 配列取得
    float* n = (float*)[mom.nBuffer contents];
    float* ux = (float*)[mom.uxBuffer contents];
    float* uy = (float*)[mom.uyBuffer contents];
    float* uz = (float*)[mom.uzBuffer contents];
    float* Pxx = (float*)[mom.PxxBuffer contents];
    float* Pxy = (float*)[mom.PxyBuffer contents];
    float* Pxz = (float*)[mom.PxzBuffer contents];
    float* Pyy = (float*)[mom.PyyBuffer contents];
    float* Pyz = (float*)[mom.PyzBuffer contents];
    float* Pzz = (float*)[mom.PzzBuffer contents];

    // グリッド情報取得
    int ngx = fldParam.ngx;
    int ngy = fldParam.ngy;
    int ngb = fldParam.ngb;
    int ng = (ngx+2*ngb+1)*(ngy+2*ngb+1);
    float dx = fldParam.dx;
    float dy = fldParam.dy;

    // minmax
    float min_n = 1e20, max_n = -1e20;
    float min_ux = 1e20, max_ux = -1e20;
    float min_uy = 1e20, max_uy = -1e20;
    float min_uz = 1e20, max_uz = -1e20;
    float min_Pxx = 1e20, max_Pxx = -1e20;
    float min_Pxy = 1e20, max_Pxy = -1e20;
    float min_Pxz = 1e20, max_Pxz = -1e20;
    float min_Pyy = 1e20, max_Pyy = -1e20;
    float min_Pyz = 1e20, max_Pyz = -1e20;
    float min_Pzz = 1e20, max_Pzz = -1e20;
    for(int i = 0; i < ng; i++){
        if (n[i]   < min_n  ) { min_n   = n[i];   } else if (n[i]   > max_n)   { max_n   = n[i]; }
        if (ux[i]  < min_ux ) { min_ux  = ux[i];  } else if (ux[i]  > max_ux)  { max_ux  = ux[i]; }
        if (uy[i]  < min_uy ) { min_uy  = uy[i];  } else if (uy[i]  > max_uy)  { max_uy  = uy[i]; }
        if (uy[i]  < min_uy ) { min_uy  = uy[i];  } else if (uy[i]  > max_uy)  { max_uy  = uy[i]; }
        if (Pxx[i] < min_Pxx) { min_Pxx = Pxx[i]; } else if (Pxx[i] > max_Pxx) { max_Pxx = Pxx[i]; }
        if (Pxy[i] < min_Pxy) { min_Pxy = Pxy[i]; } else if (Pxy[i] > max_Pxy) { max_Pxy = Pxy[i]; }
        if (Pxz[i] < min_Pxz) { min_Pxz = Pxz[i]; } else if (Pxz[i] > max_Pxz) { max_Pxz = Pxz[i]; }
        if (Pyy[i] < min_Pyy) { min_Pyy = Pyy[i]; } else if (Pyy[i] > max_Pyy) { max_Pyy = Pyy[i]; }
        if (Pyz[i] < min_Pyz) { min_Pyz = Pyz[i]; } else if (Pyz[i] > max_Pyz) { max_Pyz = Pyz[i]; }
        if (Pzz[i] < min_Pzz) { min_Pzz = Pzz[i]; } else if (Pzz[i] > max_Pzz) { max_Pzz = Pzz[i]; }
    }
    std::map<std::string, std::string> data ={
        {"n_min", fmtSci(min_n, 6)}, {"n_max", fmtSci(max_n, 6)},
        {"ux_min", fmtSci(min_ux, 6)}, {"ux_max", fmtSci(max_ux, 6)},
        {"uy_min", fmtSci(min_uy, 6)}, {"uy_max", fmtSci(max_uy, 6)},
        {"uz_min", fmtSci(min_uz, 6)}, {"uz_max", fmtSci(max_uz, 6)},
        {"Pxx_min", fmtSci(min_Pxx*0.1, 6)}, {"Pxx_max", fmtSci(max_Pxx*0.1, 6)},
        {"Pxy_min", fmtSci(min_Pxy*0.1, 6)}, {"Pxy_max", fmtSci(max_Pxy*0.1, 6)},
        {"Pxz_min", fmtSci(min_Pxz*0.1, 6)}, {"Pxz_max", fmtSci(max_Pxz*0.1, 6)},
        {"Pyy_min", fmtSci(min_Pyy*0.1, 6)}, {"Pyy_max", fmtSci(max_Pyy*0.1, 6)},
        {"Pyz_min", fmtSci(min_Pyz*0.1, 6)}, {"Pyz_max", fmtSci(max_Pyz*0.1, 6)},
        {"Pzz_min", fmtSci(min_Pzz*0.1, 6)}, {"Pzz_max", fmtSci(max_Pzz*0.1, 6)},
    };
    std::string pName_str = [pName UTF8String];
    logger.logSection("outputMoment_"+pName_str, data);

    FILE *fp = fopen(filePath, "wb");
    if (!fp) {
        NSLog(@"Error: Unable to open file %s for writing", filePath);
        return;
    }
    
    // ヘッダ情報としてメッシュ情報を出力
    fwrite(&ngx, sizeof(int), 1, fp);
    fwrite(&ngy, sizeof(int), 1, fp);
    fwrite(&ngb, sizeof(int), 1, fp);
    fwrite(&dx, sizeof(float), 1, fp);
    fwrite(&dy, sizeof(float), 1, fp);
    
    // フィールドデータを書き出す: name,type,array
    writeField(fp, [[pName stringByAppendingString:@"_n"] UTF8String], 0, n, ng, 1.0f);        // 1/cm3
    writeField(fp, [[pName stringByAppendingString:@"_ux"] UTF8String], 0, ux, ng, 1.0f);      // cm/s
    writeField(fp, [[pName stringByAppendingString:@"_uy"] UTF8String], 0, uy, ng, 1.0f);      // cm/s
    writeField(fp, [[pName stringByAppendingString:@"_uz"] UTF8String], 0, uz, ng, 1.0f);      // cm/s
    writeField(fp, [[pName stringByAppendingString:@"_Pxx"] UTF8String], 0, Pxx, ng, 0.1f);    // Pa
    writeField(fp, [[pName stringByAppendingString:@"_Pxy"] UTF8String], 0, Pxy, ng, 0.1f);    // Pa
    writeField(fp, [[pName stringByAppendingString:@"_Pxz"] UTF8String], 0, Pxz, ng, 0.1f);    // Pa
    writeField(fp, [[pName stringByAppendingString:@"_Pyy"] UTF8String], 0, Pyy, ng, 0.1f);    // Pa
    writeField(fp, [[pName stringByAppendingString:@"_Pyz"] UTF8String], 0, Pyz, ng, 0.1f);    // Pa
    writeField(fp, [[pName stringByAppendingString:@"_Pzz"] UTF8String], 0, Pzz, ng, 0.1f);    // Pa
    
    fclose(fp);
    NSLog(@"Field data successfully written to %s", filePath);
}

static void writeField(FILE* fp, const char* name, int type_id, float* array, int arrSize, float scale){
    // 変数名（固定長32文字）
    char name_buf[32] = {};
    strncpy(name_buf, name, sizeof(name_buf) - 1); // 末尾NULL保護
    fwrite(name_buf, sizeof(char), 32, fp);
    
    // タイプID（int）
    // 0: node-center value,    1: left-shifted value
    // 2: bottom-shifted value, 3: left&bottom-shifted value
    // 4: potential
    fwrite(&type_id, sizeof(int), 1, fp);
    
    // 配列本体をスケーリングして出力
    float* output = (float *)malloc(sizeof(float)*arrSize);
    for (int i = 0; i < arrSize; ++i) {
        output[i] = array[i] * scale;
    }
    fwrite(output, sizeof(float), arrSize, fp);
    // ここで必ず free する
    free(output);
}

BOOL saveProgress(int cycle, NSMutableArray* ptclArr, EMField* fld, Init* init){
    struct FlagForEquation EqFlags = init.flagForEquation;
    struct ParamForComputing compParam = init.paramForComputing;
    struct ParamForTimeIntegration timeParam = init.paramForTimeIntegration;

    // 粒子情報出力
    for (int s = 0; s < EqFlags.Particle; s++) {
        Particle *ptcl = [ptclArr objectAtIndex:s];
        SimulationParams* prm = (SimulationParams*)[[ptcl paramsBuffer] contents];
        ParticleState* p = (ParticleState*)[[ptcl particleBuffer] contents];
        NSLog(@"before: p[0].x = %e, p[0].y = %e, p[0].vx = %e",p[0].x,p[0].y,p[0].vx);
        NSLog(@"before: p[1].x = %e, p[1].y = %e, p[1].vx = %e",p[1].x,p[1].y,p[1].vx);
        NSString* fileName = [NSString stringWithFormat:@"bin/%@_ProgressData_%@_%08d.bin", timeParam.ProjectName, ptcl.pName, cycle];
        if(!saveBuffer(ptcl.particleBuffer, prm->pNum*sizeof(ParticleState), compParam.pNumMax*sizeof(ParticleState), fileName)){
            NSLog(@"[Error] output progress for %@ is failed.", ptcl.pName);
            return NO;
        }
        NSLog(@"Restart data successfully written to %@", fileName);
    }

    // 電場出力
    if (EqFlags.EMField == 1){
        uint nx = fld.ngx + 2*fld.ngb;
        uint ny = fld.ngy + 2*fld.ngb;
        NSString* fileNameEx = [NSString stringWithFormat:@"bin/%@_ProgressData_Ex_%08d.bin", timeParam.ProjectName, cycle];
        NSString* fileNameEy = [NSString stringWithFormat:@"bin/%@_ProgressData_Ey_%08d.bin", timeParam.ProjectName, cycle];
        NSString* fileNameEz = [NSString stringWithFormat:@"bin/%@_ProgressData_Ez_%08d.bin", timeParam.ProjectName, cycle];
        NSString* fileNameBx = [NSString stringWithFormat:@"bin/%@_ProgressData_Bx_%08d.bin", timeParam.ProjectName, cycle];
        NSString* fileNameBy = [NSString stringWithFormat:@"bin/%@_ProgressData_By_%08d.bin", timeParam.ProjectName, cycle];
        NSString* fileNameBz = [NSString stringWithFormat:@"bin/%@_ProgressData_Bz_%08d.bin", timeParam.ProjectName, cycle];
        if( !saveBuffer(fld.ExBuffer, (nx+2)*(ny+1)*sizeof(float), (nx+2)*(ny+1)*sizeof(float), fileNameEx) ||
            !saveBuffer(fld.EyBuffer, (nx+1)*(ny+2)*sizeof(float), (nx+1)*(ny+2)*sizeof(float), fileNameEy) ||
            !saveBuffer(fld.EzBuffer, (nx+1)*(ny+1)*sizeof(float), (nx+1)*(ny+1)*sizeof(float), fileNameEz) ||
            !saveBuffer(fld.BxBuffer, (nx+1)*(ny+2)*sizeof(float), (nx+1)*(ny+2)*sizeof(float), fileNameBx) ||
            !saveBuffer(fld.ByBuffer, (nx+2)*(ny+1)*sizeof(float), (nx+2)*(ny+1)*sizeof(float), fileNameBy) ||
            !saveBuffer(fld.BzBuffer, (nx+2)*(ny+2)*sizeof(float), (nx+2)*(ny+2)*sizeof(float), fileNameBz)){
            NSLog(@"[Error] output progress for field is failed.");
            return NO;
        }
        NSLog(@"Restart data successfully written to %@ ~ %@", fileNameEx, fileNameBz);
    }

    return YES;
}

BOOL loadProgress(int cycle, NSMutableArray* ptclArr, EMField* fld, Init* init){
    struct FlagForEquation EqFlags = init.flagForEquation;
    struct ParamForComputing compParam = init.paramForComputing;
    struct ParamForTimeIntegration timeParam = init.paramForTimeIntegration;

    // 粒子読み込み
    for (int s = 0; s < EqFlags.Particle; s++) {
        Particle *ptcl = [ptclArr objectAtIndex:s];
        SimulationParams* prm = (SimulationParams*)[[ptcl paramsBuffer] contents];
        ParticleState* p = (ParticleState*)[[ptcl particleBuffer] contents];
        NSLog(@"before: p[0].x = %e, p[0].y = %e, p[0].vx = %e",p[0].x,p[0].y,p[0].vx);
        NSLog(@"before: p[1].x = %e, p[1].y = %e, p[1].vx = %e",p[1].x,p[1].y,p[1].vx);
        NSString* fileName = [NSString stringWithFormat:@"bin/%@_ProgressData_%@_%08d.bin", timeParam.RestartName, ptcl.pName, cycle];
        uint pNum_recv;
        if(!loadBuffer(ptcl.particleBuffer, pNum_recv, compParam.pNumMax*sizeof(ParticleState), fileName)){
            NSLog(@"[Error] load progress for %@ is failed.", ptcl.pName);
            return NO;
        }
        NSLog(@"after: p[0].x = %e, p[0].y = %e, p[0].vx = %e",p[0].x,p[0].y,p[0].vx);
        NSLog(@"after: p[1].x = %e, p[1].y = %e, p[1].vx = %e",p[1].x,p[1].y,p[1].vx);
        // copy to buffer
        prm->pNum = pNum_recv/sizeof(ParticleState);
    }

    // 電場読み込み
    if (EqFlags.EMField == 1){
        uint nx = fld.ngx + 2*fld.ngb;
        uint ny = fld.ngy + 2*fld.ngb;
        uint ng_recv;
        NSString* fileNameEx = [NSString stringWithFormat:@"bin/%@_ProgressData_Ex_%08d.bin", timeParam.RestartName, cycle];
        NSString* fileNameEy = [NSString stringWithFormat:@"bin/%@_ProgressData_Ey_%08d.bin", timeParam.RestartName, cycle];
        NSString* fileNameEz = [NSString stringWithFormat:@"bin/%@_ProgressData_Ez_%08d.bin", timeParam.RestartName, cycle];
        NSString* fileNameBx = [NSString stringWithFormat:@"bin/%@_ProgressData_Bx_%08d.bin", timeParam.RestartName, cycle];
        NSString* fileNameBy = [NSString stringWithFormat:@"bin/%@_ProgressData_By_%08d.bin", timeParam.RestartName, cycle];
        NSString* fileNameBz = [NSString stringWithFormat:@"bin/%@_ProgressData_Bz_%08d.bin", timeParam.RestartName, cycle];
        if( !loadBuffer(fld.ExBuffer, ng_recv, (nx+2)*(ny+1)*sizeof(float), fileNameEx) ||
            !loadBuffer(fld.EyBuffer, ng_recv, (nx+1)*(ny+2)*sizeof(float), fileNameEy) ||
            !loadBuffer(fld.EzBuffer, ng_recv, (nx+1)*(ny+1)*sizeof(float), fileNameEz) ||
            !loadBuffer(fld.BxBuffer, ng_recv, (nx+1)*(ny+2)*sizeof(float), fileNameBx) ||
            !loadBuffer(fld.ByBuffer, ng_recv, (nx+2)*(ny+1)*sizeof(float), fileNameBy) ||
            !loadBuffer(fld.BzBuffer, ng_recv, (nx+2)*(ny+2)*sizeof(float), fileNameBz)){
            NSLog(@"[Error] load progress for field is failed.");
            return NO;
        }
    }

    return YES;
}

//------------------------------------------------------------------------------
// Metal バッファの内容をバイナリファイルに書き出すメソッド
//
// @param buffer      書き出したい MTLBuffer
// @param length      書き出す要素数
// @param lengthMax   確保中のバッファサイズ
// @param outPath     保存先のファイルパス (NSString*)
// @return YES: 成功, NO: 失敗
//------------------------------------------------------------------------------
BOOL saveBuffer(id<MTLBuffer> buffer, uint length, uint lengthMax, NSString* outPath){
    // 0) サイズチェック
    if(length > lengthMax){
        NSLog(@"[Error] 書き出しサイズが確保中の配列サイズを超えています。");
        return NO;
    }

    // 1) MTLBuffer の内容を CPU アクセス可能な生ポインタとして取得
    void *rawPtr = [buffer contents];
    if (!rawPtr) {
        NSLog(@"[Error] buffer.contents が NULL です。");
        return NO;
    }

    // 2) ファイルをバイナリ書き込みモードでオープン
    const char *cPath = [outPath fileSystemRepresentation];
    FILE *fp = fopen(cPath, "wb");
    if (!fp) {
        NSLog(@"[Error] ファイルオープンに失敗しました: %s", cPath);
        return NO;
    }

    // 3) 生ポインタから length 要素を書き込む
    size_t written = fwrite(rawPtr, 1, length, fp);
    fclose(fp);

    if (written != length) {
        NSLog(@"[Error] fwrite の要素数が一致しません (written: %zu, expected: %lu)", written, (unsigned long)length);
        return NO;
    }

    return YES;
}

//------------------------------------------------------------------------------
// Metal バッファの内容をバイナリファイルから読み込むメソッド
//
// @param buffer      読み込んだ内容を格納する MTLBuffer
// @param length      読み込む要素数
// @param lengthMax   確保中のバッファサイズ
// @param inPath      読み込み元のファイルパス (NSString*)
// @return YES: 成功, NO: 失敗
//------------------------------------------------------------------------------
BOOL loadBuffer(id<MTLBuffer> buffer, uint& length, uint lengthMax, NSString* inPath){
    // 1) ファイルをバイナリ読み込みモードでオープン
    const char *cPath = [inPath fileSystemRepresentation];
    FILE *fp = fopen(cPath, "rb");
    if (!fp) {
        NSLog(@"[Error] ファイルオープンに失敗しました: %s", cPath);
        return NO;
    }

    // 2) ファイルサイズを取得して length に格納
    fseek(fp, 0, SEEK_END);
    long fileSize = ftell(fp);
    rewind(fp);
    if (fileSize < 0) {
        NSLog(@"[Error] ftell でファイルサイズ取得に失敗");
        fclose(fp);
        return NO;
    }
    length = fileSize;

    // 3) サイズチェック
    if(length > lengthMax){
        NSLog(@"[Error] 読み込みサイズが確保中の配列サイズを超えています。");
        fclose(fp);
        return NO;
    }
    
    // 4) MTLBuffer の内容を CPU アクセス可能な生ポインタとして取得
    void *rawPtr = [buffer contents];
    if (!rawPtr) {
        NSLog(@"[Error] buffer.contents が NULL です。");
        fclose(fp);
        return NO;
    }

    // 4) ファイルから rawPtr へ一気に読み込む
    size_t readBytes = fread(rawPtr, 1, fileSize, fp);
    fclose(fp);


    return YES;
}