#include <Python.h>

int main() {
    // Initialize the Python Interpreter
    Py_Initialize();

    // Add the current directory to the Python path
    PyRun_SimpleString("import sys; sys.path.append('/Users/mofefagade/CLionProjects/Call_1')");

    // Import the Python module 'mytest'
    PyRun_SimpleString("import mytest");

    // Call a function from 'mytest' and print the result
    PyRun_SimpleString("print(mytest.myabs(2.0))");

    // Finalize the Python Interpreter
    Py_Finalize();

    return 0;
}
