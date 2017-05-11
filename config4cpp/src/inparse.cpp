#include "inparse.h"



namespace rockseis {

// constructor
Inparse::Inparse()
{
    cfg =  config4cpp::Configuration::create();
    scope = (char *) calloc(1,1);
}

bool Inparse::parse(std::string configfile)
{
    configFile = configfile;
    try {
        cfg->parse(configFile.c_str());
    } catch(const config4cpp::ConfigurationException & ex) {
        std::cerr << ex.c_str() << std::endl;
        return INPARSE_ERR;
    }
    return INPARSE_OK;
}

bool Inparse::getPar(std::string par, int *var)
{
    try {
        *var = cfg->lookupInt(scope, par.c_str());
    } catch(const config4cpp::ConfigurationException & ex) {
        std::cerr << ex.c_str() << std::endl;
        return INPARSE_ERR;
    }
    return INPARSE_OK;
}

bool Inparse::getPar(std::string par, float *var)
{
    try {
        *var = cfg->lookupFloat(scope, par.c_str());
    } catch(const config4cpp::ConfigurationException & ex) {
        std::cerr << ex.c_str() << std::endl;
        return INPARSE_ERR;
    }
    return INPARSE_OK;
}

bool Inparse::getPar(std::string par, bool *var)
{
    try {
        *var = cfg->lookupBoolean(scope, par.c_str());
    } catch(const config4cpp::ConfigurationException & ex) {
        std::cerr << ex.c_str() << std::endl;
        return INPARSE_ERR;
    }
    return INPARSE_OK;
}

bool Inparse::getPar(std::string par, std::string *var)
{
    try {
        *var = cfg->lookupString(scope, par.c_str());
    } catch(const config4cpp::ConfigurationException & ex) {
        std::cerr << ex.c_str() << std::endl;
        return INPARSE_ERR;
    }
    return INPARSE_OK;
}

// Destructor
Inparse::~Inparse() {
	// Destroy cfg
	cfg->destroy();
}

}
