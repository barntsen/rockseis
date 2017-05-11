#ifndef INPARSE_H
#define INPARSE_H

// Include statements
#include <iostream>
#include <vector>
#include <string>
#include <memory>
#include <config4cpp/Configuration.h>

#define INPARSE_OK 1
#define INPARSE_ERR 0

namespace rockseis {

// =============== INPARSE CLASS =============== //
/** The Inparse class
 *
 */
class Inparse {
public:
    Inparse(); ///<Constructor
    virtual ~Inparse();	///< Destructor
    bool parse(std::string configfile);
    bool getPar(std::string par, int *var); ///< Get int
    bool getPar(std::string par, float *var); ///< Get float
    bool getPar(std::string par, bool *var); ///< Get boolean
    bool getPar(std::string par, std::string *var); ///< Get string

private:
    config4cpp::Configuration *  cfg;
    const char *     scope;
    std::string      configFile;
};
}
#endif //INPARSE_H
