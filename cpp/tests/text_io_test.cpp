#include "util/text_stream.hpp"
#include <fstream>
#include <iterator>
#include <iostream>
#include <string>
int main() {
 const std::string p="tmp/preserved-result.txt";
 { std::ofstream f(p,std::ios::binary);f<<"existing result\n"; }
 TextWriter::assertWritable(p);
 std::ifstream f(p,std::ios::binary);
 std::string s((std::istreambuf_iterator<char>(f)),std::istreambuf_iterator<char>());
 if(s!="existing result\n")return 1;
 f.close();std::remove(p.c_str());
 {TextWriter w(p);w.write("one\ntwo\n");w.close();}
 std::ifstream g(p,std::ios::binary);s.assign(std::istreambuf_iterator<char>(g),{});
 g.close();std::remove(p.c_str());
 if(s!="one\ntwo\n")return 2;
 std::cout<<"Text I/O: existing-output preservation and canonical LF passed.\n";
}

