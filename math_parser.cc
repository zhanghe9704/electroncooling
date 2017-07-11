#include "math_parser.h"

#include <cassert>
#include <iostream>
#include <string>

#include "constants.h"
#include "muParserDLL.h"

#define PARSER_MAXVARS		500

#ifndef _UNICODE
    #define _T(x) x
    #define myprintf printf
    #define mystrlen strlen
    #define myfgets fgets
    #define mystrcmp strcmp
#else
#define _T(x) L ##x
    #define myprintf wprintf
    #define mystrlen wcslen
    #define myfgets fgetws
    #define mystrcmp wcscmp
#endif

//---------------------------------------------------------------------------
// Factory function for creating new parser variables
// This could as well be a function performing database queries.
muFloat_t* AddVariable(const muChar_t* a_szName, void *pUserData)
{
    static muFloat_t afValBuf[PARSER_MAXVARS];  // I don't want dynamic allocation here
    static int iVal = 0;                     // so i used this buffer

//    myprintf(_T("Generating new variable \"%s\" (slots left: %d; context pointer: 0x%x)\n"), a_szName, PARSER_MAXVARS - iVal, (int)pUserData);
//    myprintf(_T("Generating new variable \"%s\" (slots left: %d \n"), a_szName, PARSER_MAXVARS - iVal);

    afValBuf[iVal] = 0;
    if (iVal >= PARSER_MAXVARS - 1)
    {
        myprintf(_T("Variable buffer overflow."));
        return NULL;
    }

    return &afValBuf[iVal++];
}

//---------------------------------------------------------------------------
void ListVar(muParserHandle_t a_hParser)
{
    int iNumVar = mupGetVarNum(a_hParser);
    int i = 0;

    if (iNumVar == 0)
    {
        myprintf(_T("No variables defined\n"));
        return;
    }

    myprintf(_T("\nExpression variables:\n"));
    myprintf(_T("---------------------\n"));
    myprintf(_T("Number: %d\n"), iNumVar);

    for (i = 0; i < iNumVar; ++i)
    {
        const muChar_t* szName = 0;
        muFloat_t* pVar = 0;

        mupGetVar(a_hParser, i, &szName, &pVar);
        myprintf(_T("Name: %s    Address: [0x%x]\n"), szName, (long long)pVar);
    }
}

//---------------------------------------------------------------------------
void ListExprVar(muParserHandle_t a_hParser)
{
    muInt_t iNumVar = mupGetExprVarNum(a_hParser),
        i = 0;

    if (iNumVar == 0)
    {
        myprintf(_T("Expression dos not contain variables\n"));
        return;
    }

    myprintf(_T("\nExpression variables:\n"));
    myprintf(_T("---------------------\n"));
    myprintf(_T("Expression: %s\n"), mupGetExpr(a_hParser));
    myprintf(_T("Number: %d\n"), iNumVar);

    for (i = 0; i < iNumVar; ++i)
    {
        const muChar_t* szName = 0;
        muFloat_t* pVar = 0;

        mupGetExprVar(a_hParser, i, &szName, &pVar);
        myprintf(_T("Name: %s   Address: [0x%x]\n"), szName, (long long)pVar);
    }
}

//---------------------------------------------------------------------------
void ListConst(muParserHandle_t a_hParser)
{
    muInt_t iNumVar = mupGetConstNum(a_hParser),
        i = 0;

    if (iNumVar == 0)
    {
        myprintf(_T("No constants defined\n"));
        return;
    }

    myprintf(_T("\nParser constants:\n"));
    myprintf(_T("---------------------\n"));
    myprintf(_T("Number: %d"), iNumVar);

    for (i = 0; i < iNumVar; ++i)
    {
        const muChar_t* szName = 0;
        muFloat_t fVal = 0;

        mupGetConst(a_hParser, i, &szName, &fVal);
        myprintf(_T("  %s = %f\n"), szName, fVal);
    }
}



//---------------------------------------------------------------------------
// Callback function for parser errors
void OnError(muParserHandle_t hParser)
{
    myprintf(_T("\nError:\n"));
    myprintf(_T("------\n"));
    myprintf(_T("Message:  \"%s\"\n"), mupGetErrorMsg(hParser));
    myprintf(_T("Token:    \"%s\"\n"), mupGetErrorToken(hParser));
    myprintf(_T("Position: %d\n"), mupGetErrorPos(hParser));
    myprintf(_T("Errc:     %d\n"), mupGetErrorCode(hParser));
    assert(false&&"ERROR IN MATH PARSING!");
}

void AddConst(muParserHandle_t hParser) {
    mupDefineConst(hParser, "k_c", k_c);
    mupDefineConst(hParser, "k_e", k_e);
    mupDefineConst(hParser, "k_pi", k_pi);
    mupDefineConst(hParser, "k_u", k_u);
    mupDefineConst(hParser, "k_me", k_me);
}

//Initialize the math parser environment
void initialize_parser(muParserHandle_t &math_parser) {
    math_parser = mupCreate(muBASETYPE_FLOAT);              // initialize the parser
    mupSetErrorHandler(math_parser, OnError);
    mupSetVarFactory(math_parser, AddVariable, NULL);       // Set a variable factory
    AddConst(math_parser);                                  // Add constants
}



int use_parser(int argc, char* argv[])
{
    muChar_t szLine[100];
    muFloat_t fVal = 0;
    muParserHandle_t hParser;
    hParser = mupCreate(muBASETYPE_FLOAT);              // initialize the parser
    mupSetErrorHandler(hParser, OnError);

    // Set a variable factory
    mupSetVarFactory(hParser, AddVariable, NULL);

    AddConst(hParser);


    ListConst(hParser);

    std::string str = "x = 10.6e-1";
    mupSetExpr(hParser, str.c_str());

    fVal = mupEval(hParser);

    ListVar(hParser);
    ListExprVar(hParser);

    return 0;
}
