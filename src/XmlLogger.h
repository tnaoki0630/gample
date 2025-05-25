#pragma once

#include <fstream>
#include <sstream>
#include <map>
#include <string>

class XmlLogger {
public:
    // コンストラクタ：出力ファイル名を受け取り、XMLヘッダを書き出す
    explicit XmlLogger(const std::string& filename) {
        ofs_.open(filename);
        if (!ofs_) throw std::runtime_error("Cannot open log file: " + filename);
        ofs_ << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
             << "<SimulationLog>\n";
    }

    // サイクル開始
    void logCycleStart(int cycle, double time){
        ofs_ << " <Cycle ID=\"" << cycle << "\">\n"
             << "  time = " << time <<"\n";
    }
    // ログ出力
    void logSection(const std::string& name, const std::map<std::string, std::string>& data) {
        ofs_ << "  <Section Name=\"" << name << "\">\n";
        for (const auto& kv : data){
            ofs_ << "   <" << kv.first << ">"
            << kv.second
            << "</" << kv.first << ">\n";
        }
        ofs_ << "  </Section>\n";
    }
    // サイクル終了
    void logCycleEnd(){
        ofs_ << " </Cycle>\n";
    }

    // コメント出力
    void logComment(const std::string& comment){
        ofs_ << "<![CDATA[\n" << comment << "]]>\n";
    }

    // デストラクタ：ルート要素を閉じてファイルを自動クローズ
    ~XmlLogger() {
        ofs_ << "</SimulationLog>\n";
    }

    // コピー禁止・ムーブ禁止(書き込み競合回避)
    XmlLogger(const XmlLogger&) = delete;
    XmlLogger& operator=(const XmlLogger&) = delete;

private:
    std::ofstream ofs_;
};

// フォーマット変換
static std::string fmtSci(float v, int prec) {
    std::ostringstream oss;
    oss << std::scientific << std::setprecision(prec) << v;
    return oss.str();
}