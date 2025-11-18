#pragma once

#import <fstream>
#import <sstream>
#import <map>
#import <string>

class XmlLogger {
public:
    // コンストラクタ：出力ファイル名を受け取り、ヘッダを書き込む
    explicit XmlLogger(const std::string& logName, const std::string& CSVheader) {
        // XML log
        _ofs.open(logName+"_log.xml");
        if (!_ofs) throw std::runtime_error("Cannot open log file: " + logName);
        _ofs << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
             << "<SimulationLog>\n";
        // time-series csv
        _ofsCSV.open(logName+"_timeseries.csv");
        if (!_ofsCSV) throw std::runtime_error("Cannot open csv file: " + logName);
        _ofsCSV << CSVheader << "\n";
        // flag
        _outputFlag = false;
    }

    // 出力スイッチ
    void swichLog(BOOL swich){
        _outputFlag = swich;
    }

    // サイクル開始
    void logCycleStart(int cycle, double time){
        if (_outputFlag) {
            // XML
            _ofs << " <Cycle ID=\"" << cycle << "\">\n"
                 << "  time = " << time <<"\n";
            // CSV
            _ofsCSV << cycle << "," << time;
        }
    }
    // ログ出力
    void logSection(const std::string& name, const std::map<std::string, std::string>& data) {
        if (_outputFlag) {
            _ofs << "  <Section Name=\"" << name << "\">\n";
            for (const auto& kv : data){
                _ofs << "   <" << kv.first << ">"
                << kv.second
                << "</" << kv.first << ">\n";
            }
            _ofs << "  </Section>\n";
        }
    }
    // csv出力
    void csvOutput(const std::string& data) {
        if (_outputFlag) {
            _ofsCSV << "," << data;
        }
    }
    // サイクル終了
    void logCycleEnd(){
        if (_outputFlag) {
            // XML
            _ofs << " </Cycle>\n";
            // CSV
            _ofsCSV << "\n";
        }
    }

    // コメント出力
    void logComment(const std::string& comment){
        _ofs << "<![CDATA[\n" << comment << "]]>\n";
    }

    // デストラクタ：ルート要素を閉じてファイルを自動クローズ
    ~XmlLogger() {
        _ofs << "</SimulationLog>\n";
    }

    // コピー禁止・ムーブ禁止(書き込み競合回避)
    XmlLogger(const XmlLogger&) = delete;
    XmlLogger& operator=(const XmlLogger&) = delete;

private:
    std::ofstream _ofs;
    std::ofstream _ofsCSV;
    BOOL _outputFlag;
};

// フォーマット変換
static std::string fmtSci(float v, int prec) {
    std::ostringstream oss;
    oss << std::scientific << std::setprecision(prec) << v;
    return oss.str();
}
