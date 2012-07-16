#ifndef CONFIG_H
#define CONFIG_H 
#include"baseconfig.h"

class CConfig : public CBaseConfig{
	public:
	CConfig(string filename){
		define_parameters();
		parse(filename);
		}
	void define_parameters(){
		add_param<size_t>("rmax",10);
		}
	};

#endif /* CONFIG_H */
