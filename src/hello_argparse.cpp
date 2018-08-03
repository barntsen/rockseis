#include <iostream>
#include <args.hxx>
int main(int argc, char **argv)
{
    args::ArgumentParser parser("This is a test program.", "This goes after the options.");
    parser.LongPrefix("-");
    parser.LongSeparator("=");
    args::HelpFlag help(parser, "help", "Display this help menu", {"h", "help"});
    args::ValueFlag<int> integer(parser, "integer", "The integer flag", {"i"});
    args::ValueFlag<std::string> input(parser, "BLOCK SIZE", "Block size", {"if"});
    args::ValueFlag<std::string> output(parser, "BLOCK SIZE", "Block size", {"of"});
    try
    {
        parser.ParseCLI(argc, argv);
    }
    catch (args::Help)
    {
        std::cout << parser;
        return 0;
    }
    catch (args::ParseError e)
    {
        std::cerr << e.what() << std::endl;
        std::cerr << parser;
        return 1;
    }
    catch (args::ValidationError e)
    {
        std::cerr << e.what() << std::endl;
        std::cerr << parser;
        return 1;
    }
    if (integer) { std::cout << "i: " << args::get(integer) << std::endl; }
    if (input) { std::cout << "if = " << args::get(input) << std::endl; }
    if (output) { std::cout << "of = " << args::get(output) << std::endl; }
    float i = args::get(integer);
    std::cerr << "i: " << i << std::endl;
    return 0;
}
