// MemoryUtils.h
#import <mach/mach.h>
#import <Foundation/Foundation.h>
#import <chrono>
#import <string>
#import <map>

namespace MemoryUtils {
    /// 現在プロセスが使っている物理メモリ（バイト）を返す
    static size_t currentPhysicalFootprint() {
        task_vm_info_data_t vmInfo;
        mach_msg_type_number_t count = TASK_VM_INFO_COUNT;
        kern_return_t err = task_info(mach_task_self(),
                                      TASK_VM_INFO,
                                      reinterpret_cast<task_info_t>(&vmInfo),
                                      &count);
        if (err != KERN_SUCCESS) {
            return 0; // 取得失敗時は 0 を返す
        }
        return static_cast<size_t>(vmInfo.phys_footprint);
    }

    /// 現在プロセスのRSS (Resident Set Size) を返す（参考情報）
    static size_t currentRSS() {
        task_basic_info_data_t info;
        mach_msg_type_number_t size = TASK_BASIC_INFO_COUNT;
        kern_return_t kr = task_info(mach_task_self(),
                                     TASK_BASIC_INFO,
                                     reinterpret_cast<task_info_t>(&info),
                                     &size);
        if (kr != KERN_SUCCESS) {
            return 0;
        }
        return static_cast<size_t>(info.resident_size);
    }
}

// tyoe of name: std::string
// tyoe of data: std::map<std::string,std::string> 
#define MEASURE(name, code, data)                                             \
    do {                                                                       \
        auto _start = std::chrono::high_resolution_clock::now();              \
        code;                                                                 \
        auto _end = std::chrono::high_resolution_clock::now();                \
        auto _time =                                                          \
            std::chrono::duration_cast<std::chrono::microseconds>(_end - _start).count(); \
        (data)[(name)] = std::to_string(_time);                                \
    } while (0)

