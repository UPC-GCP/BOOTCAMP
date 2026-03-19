
#include <iostream>
#include <numbers>

#define exprtk_enable_all_features
#define exprtk_disable_string_capabilities
#include "exprtk.hpp"



// Lo estoy haciendo como una función por ahora para tenerlo implementado y hacer funcionar el parser. Después creo que sería buena idea tener class separado y que MainSolver.cpp se encargue de cargarle las variables para no tener que volver a crear el symbol table cada vez

inline double calcFormulaExpression(std::string exprFormula, double exprVal){

    // Variables
    exprtk::symbol_table<double> symbol_table;
    exprtk::expression<double> expression;
    exprtk::parser<double> parser;
    double TVar = exprVal, piVar = M_PI;

    // Parsing
    symbol_table.add_variable("t", TVar);
    symbol_table.add_variable("pi", piVar);
    expression.register_symbol_table(symbol_table);
    parser.compile(exprFormula, expression);

    return expression.value();

}
