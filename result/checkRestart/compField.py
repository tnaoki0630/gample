import numpy as np
import matplotlib.pyplot as plt

def getFields(filename):
    with open(filename, 'rb') as f:
        # ヘッダ情報の読み込み
        ngx = np.frombuffer(f.read(4), dtype=np.int32)[0]
        ngy = np.frombuffer(f.read(4), dtype=np.int32)[0]
        ngb = np.frombuffer(f.read(4), dtype=np.int32)[0]
        dx = np.frombuffer(f.read(4), dtype=np.float32)[0]
        dy = np.frombuffer(f.read(4), dtype=np.float32)[0]

        print(f"Mesh info: ngx={ngx}, ngy={ngy}, ngb={ngb}, dx={dx}, dy={dy}")

        fields = []
        while True:
            # --- 1. 変数名（固定長32バイト）
            name_bytes = f.read(32)
            if not name_bytes:
                break  # EOF
            # 終端のNULL文字を削除して文字列として格
            name = name_bytes.decode("utf-8").rstrip('\x00')

            # --- 2. タイプID（int）
            type_id = int.from_bytes(f.read(4), byteorder="little")

            # 配列サイズの決定
            nx = ngx+2*ngb
            ny = ngy+2*ngb
            if type_id == 4:
                nx +=1
                ny +=1
            
            # --- 3. 配列
            arrSize = (nx+1)*(ny+1)
            arr = np.frombuffer(f.read(arrSize * 4), dtype=np.float32).reshape((ny+1, nx+1))

            # タプルに格納
            fields.append((name, type_id, arr))
            print(name, type_id, arrSize, arr)
    return fields

if __name__ == '__main__':

    cycle = 600
    filename = f"bin/proj1_EMField_{cycle:08}.bin"
    filename = f"bin/proj1_Moments_electron_{cycle:08}.bin"
    filename = f"bin/proj1_Moments_ion_Xe1_{cycle:08}.bin"
    fields = getFields(filename)
    filename = f"bin/proj1_restart_EMField_{cycle:08}.bin"
    filename = f"bin/proj1_restart_Moments_electron_{cycle:08}.bin"
    filename = f"bin/proj1_restart_Moments_ion_Xe1_{cycle:08}.bin"
    fields_restart = getFields(filename)
    
    for (name, type_id, arr), (name_r, type_id_r, arr_r) in zip(fields, fields_restart):
        diff = abs(arr-arr_r)/abs(arr)
        print(f"[org]     name: {name}, arr.min = {arr.min()}, arr.max = {arr.max()}")
        print(f"[restart] name: {name_r}, arr.min = {arr_r.min()}, arr.max = {arr_r.max()}")
        print(f"[diff]    name: {name}, diff.min = {diff.min()}, diff.max = {diff.max()}")
    
