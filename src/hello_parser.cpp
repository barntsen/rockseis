#include <config4cpp/Configuration.h>
#include <iostream>
using namespace config4cpp;
using namespace std;

int main(int argc, char ** argv)
{
    Configuration *  cfg = Configuration::create();
    const char *     scope = "";
    const char *     configFile = "configfile.cfg";
    const char *     url;
    const char *     file;
    bool             true_false;
    float            somenumber;

    try {
        cfg->parse(configFile);
        url        = cfg->lookupString(scope, "url");
        file       = cfg->lookupString(scope, "file");
        true_false = cfg->lookupBoolean(scope, "true_false");
        somenumber = cfg->lookupFloat(scope, "somenumber");
    } catch(const ConfigurationException & ex) {
        cerr << ex.c_str() << endl;
        cfg->destroy();
        return 1;
    }
    cout << "url=" << url << "; file=" << file
         << "; true_false=" << true_false
         << "; somenumber=" << somenumber
         << endl;
    cfg->destroy();
    return 0;
}

