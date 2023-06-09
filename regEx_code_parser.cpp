#include <iostream>
#include <regex>
#include <fstream>
#include <string>

void identifyFunctions(const std::string& code) {
    std::regex functionRegex(R"(\b([a-zA-Z_][a-zA-Z0-9_]*)\s*\([^)]*\)\s*\{)");

    std::sregex_iterator regexIterator(code.begin(), code.end(), functionRegex);
    std::sregex_iterator regexIteratorEnd;

    while (regexIterator != regexIteratorEnd) {
        std::smatch match = *regexIterator;
        std::string functionName = match.str(1);
        std::cout << "Function: " << functionName << std::endl;

        ++regexIterator;
    }
}

int main() {
    std::ifstream file("path/to/your/cpp/file.cpp");
    if (!file) {
        std::cerr << "Failed to open file." << std::endl;
        return 1;
    }

    std::string code((std::istreambuf_iterator<char>(file)),
                     std::istreambuf_iterator<char>());

    identifyFunctions(code);

    return 0;
}

