#import "FileIO.h"
BOOL saveProgress(int cycle, NSMutableArray* ptclArr, EMField* fld, int current, Init* init){
    struct FlagForEquation EqFlags = init.flagForEquation;
    struct ParamForComputing compParam = init.paramForComputing;

    // 粒子情報出力
    for (int s = 0; s < EqFlags.Particle; s++) {
        Particle *ptcl = [ptclArr objectAtIndex:s];
        SimulationParams* prm = (SimulationParams*)[[ptcl paramsBuffer] contents];
        NSString* fileName = [NSString stringWithFormat:@"bin/restart_%@_%08d.bin", ptcl.pName, cycle];
        if(!saveBuffer(ptcl.particleBuffer, prm->pNum, compParam.pNumMax, fileName)){
            NSLog(@"[Error] output progress for %@ is failed.", ptcl.pName);
            return NO;
        }
        NSLog(@"Restart data successfully written to %@", fileName);
    }

    // 電場出力
    if (EqFlags.EMField == 1){
        uint ng = (fld.ngx + 2*fld.ngb + 1)*(fld.ngy + 2*fld.ngb + 1);
        NSString* fileNameEx = [NSString stringWithFormat:@"bin/restart_Ex_%08d.bin", cycle];
        NSString* fileNameEy = [NSString stringWithFormat:@"bin/restart_Ey_%08d.bin", cycle];
        NSString* fileNameEz = [NSString stringWithFormat:@"bin/restart_Ez_%08d.bin", cycle];
        NSString* fileNameBx = [NSString stringWithFormat:@"bin/restart_Bx_%08d.bin", cycle];
        NSString* fileNameBy = [NSString stringWithFormat:@"bin/restart_By_%08d.bin", cycle];
        NSString* fileNameBz = [NSString stringWithFormat:@"bin/restart_Bz_%08d.bin", cycle];
        if( !saveBuffer(fld.ExBuffer, ng, ng, fileNameEx) ||
            !saveBuffer(fld.EyBuffer, ng, ng, fileNameEy) ||
            !saveBuffer(fld.EzBuffer, ng, ng, fileNameEz) ||
            !saveBuffer(fld.BxBuffer, ng, ng, fileNameBx) ||
            !saveBuffer(fld.ByBuffer, ng, ng, fileNameBy) ||
            !saveBuffer(fld.BzBuffer, ng, ng, fileNameBz)){
            NSLog(@"[Error] output progress for field is failed.");
            return NO;
        }
        NSLog(@"Restart data successfully written to %@ ~ %@", fileNameEx, fileNameBz);
    }

    return YES;
}

BOOL loadProgress(int cycle, NSMutableArray* ptclArr, EMField* fld, int& current, Init* init){
    struct FlagForEquation EqFlags = init.flagForEquation;
    struct ParamForComputing compParam = init.paramForComputing;

    // 粒子読み込み
    for (int s = 0; s < EqFlags.Particle; s++) {
        Particle *ptcl = [ptclArr objectAtIndex:s];
        SimulationParams* prm = (SimulationParams*)[[ptcl paramsBuffer] contents];
        NSString* fileName = [NSString stringWithFormat:@"bin/restart_%@_%08d.bin", ptcl.pName, cycle];
        uint pNum_recv;
        if(!loadBuffer(ptcl.particleBuffer, pNum_recv, compParam.pNumMax, fileName)){
            NSLog(@"[Error] load progress for %@ is failed.", ptcl.pName);
            return NO;
        }
        // copy to buffer
        prm->pNum = pNum_recv;
    }

    // 電場読み込み
    if (EqFlags.EMField == 1){
        uint ng = (fld.ngx + 2*fld.ngb + 1)*(fld.ngy + 2*fld.ngb + 1);
        uint ng_recv;
        NSString* fileNameEx = [NSString stringWithFormat:@"bin/restart_Ex_%08d.bin", cycle];
        NSString* fileNameEy = [NSString stringWithFormat:@"bin/restart_Ey_%08d.bin", cycle];
        NSString* fileNameEz = [NSString stringWithFormat:@"bin/restart_Ez_%08d.bin", cycle];
        NSString* fileNameBx = [NSString stringWithFormat:@"bin/restart_Bx_%08d.bin", cycle];
        NSString* fileNameBy = [NSString stringWithFormat:@"bin/restart_By_%08d.bin", cycle];
        NSString* fileNameBz = [NSString stringWithFormat:@"bin/restart_Bz_%08d.bin", cycle];
        if( !loadBuffer(fld.ExBuffer, ng_recv, ng, fileNameEx) ||
            !loadBuffer(fld.EyBuffer, ng_recv, ng, fileNameEy) ||
            !loadBuffer(fld.EzBuffer, ng_recv, ng, fileNameEz) ||
            !loadBuffer(fld.BxBuffer, ng_recv, ng, fileNameBx) ||
            !loadBuffer(fld.ByBuffer, ng_recv, ng, fileNameBy) ||
            !loadBuffer(fld.BzBuffer, ng_recv, ng, fileNameBz)){
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
    size_t written = fwrite(rawPtr, 1, length*sizeof(float), fp);
    fclose(fp);

    if (written != length*sizeof(float)) {
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
    length = fileSize/sizeof(float);

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