/*

Copyright (c) 2014 Jarryd Beck

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

*/

#include <iostream>

#include "cxxopts.hpp"

int main(int argc, char* argv[])
{
  try
  {
    cxxopts::Options options(argv[0], " - example command line options");
    options
      .positional_help("[optional args]")
      .show_positional_help();

    bool apple = false;

    options
      .allow_unrecognised_options()
      .add_options()
      ("a,apple", "an apple", cxxopts::value<bool>(apple))
      ("b,bob", "Bob")
      ("t,true", "True", cxxopts::value<bool>()->default_value("true"))
      ("f, file", "File", cxxopts::value<std::vector<std::string>>(), "FILE")
      ("i,input", "Input", cxxopts::value<std::string>())
      ("o,output", "Output file", cxxopts::value<std::string>()
          ->default_value("a.out")->implicit_value("b.def"), "BIN")
      ("positional",
        "Positional arguments: these are the arguments that are entered "
        "without an option", cxxopts::value<std::vector<std::string>>())
      ("long-description",
        "thisisareallylongwordthattakesupthewholelineandcannotbebrokenataspace")
      ("help", "Print help")
      ("int", "An integer", cxxopts::value<int>(), "N")
      ("float", "A floating point number", cxxopts::value<float>())
      ("option_that_is_too_long_for_the_help", "A very long option")
    #ifdef CXXOPTS_USE_UNICODE
      ("unicode", u8"A help option with non-ascii: à. Here the size of the"
        " string should be correct")
    #endif
    ;

    options.add_options("Group")
      ("c,compile", "compile")
      ("d,drop", "drop", cxxopts::value<std::vector<std::string>>());

    options.parse_positional({"input", "output", "positional"});

    auto result = options.parse(argc, argv);

    if (result.count("help"))
    {
      std::cerr << options.help({"", "Group"}) << std::endl;
      exit(0);
    }

    if (apple)
    {
      std::cerr << "Saw option ‘a’ " << result.count("a") << " times " <<
        std::endl;
    }

    if (result.count("b"))
    {
      std::cerr << "Saw option ‘b’" << std::endl;
    }

    if (result.count("f"))
    {
      auto& ff = result["f"].as<std::vector<std::string>>();
      std::cerr << "Files" << std::endl;
      for (const auto& f : ff)
      {
        std::cerr << f << std::endl;
      }
    }

    if (result.count("input"))
    {
      std::cerr << "Input = " << result["input"].as<std::string>()
        << std::endl;
    }

    if (result.count("output"))
    {
      std::cerr << "Output = " << result["output"].as<std::string>()
        << std::endl;
    }

    if (result.count("positional"))
    {
      std::cerr << "Positional = {";
      auto& v = result["positional"].as<std::vector<std::string>>();
      for (const auto& s : v) {
        std::cerr << s << ", ";
      }
      std::cerr << "}" << std::endl;
    }

    if (result.count("int"))
    {
      std::cerr << "int = " << result["int"].as<int>() << std::endl;
    }

    if (result.count("float"))
    {
      std::cerr << "float = " << result["float"].as<float>() << std::endl;
    }

    std::cerr << "Arguments remain = " << argc << std::endl;

  } catch (const cxxopts::OptionException& e)
  {
    std::cerr << "error parsing options: " << e.what() << std::endl;
    exit(1);
  }

  return 0;
}
