#ifndef MATH_PARSER_H_INCLUDED
#define MATH_PARSER_H_INCLUDED

#include "muParserDLL.h"

void ListVar(muParserHandle_t a_hParser);
void ListExprVar(muParserHandle_t a_hParser);
void ListConst(muParserHandle_t a_hParser);
void initialize_parser(muParserHandle_t &math_parser);



#endif // MATH_PARSER_H_INCLUDED
