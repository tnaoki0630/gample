// main.mm
#include <iostream>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
using namespace boost::property_tree;
using namespace boost::property_tree::xml_parser;

int main(int argc, const char * argv[]) {
    // XMLファイルを読み込む
	ptree pt;
	read_xml("CrossSection_hayashi_e_Xe.xml", pt, trim_whitespace);

    // 子要素を取得
	if (auto opt = pt.get_child_optional("root.sample0"))
	{	
		for (auto& child : opt.value())
		{
            // 要素名 取得
			std::cout << child.first.c_str() << std::endl;	

			// 属性color 取得
			if (auto opt = child.second.get_optional<std::string>("<xmlattr>.color"))			
			{
				std::cout << "\t" << opt.value() << std::endl;
			}
			else
			{
				// 指定した属性がなかったとき
				std::cout << "\t" << "属性なし" << std::endl;
			}

			// 属性id 取得
			if (auto opt = child.second.get_optional<std::string>("<xmlattr>.id"))
			{
				std::cout << "\t" << opt.value() << std::endl;
			}
            // elseは省略   

			// 要素の値
			if (auto opt = child.second.get_value_optional<std::string>())
			{
				std::cout << "\t" << opt.value() << std::endl;
			}
			
			// さらに子要素を取得			
			if (auto opt = child.second.get_child_optional(""))
			{
				opt.value().erase("<xmlattr>");		// 属性情報を削除
				for (auto& item : opt.value())
				{
					std::cout << "\t" << item.second.get_value<std::string>() << std::endl;
				}
			}
		}
	}
	else
	{
		std::cout << "<root>.<sample0>が存在しない." << std::endl;
	}

    return 0;
}
