#include <Python.h>
#include <stdio.h>

int main() {
    // Initialize the Python Interpreter
    Py_Initialize();

    // Print initial message
    printf("Running the program...\n");

    // Add current directory to Python path
    PyRun_SimpleString("import sys; sys.path.append(\".\")");

    // Set the name of the Python module
    PyObject *pName = PyUnicode_DecodeFSDefault("numpy_example");

    // Import the module
    PyObject *pModule = PyImport_Import(pName);
    Py_DECREF(pName);

    // Check if the module was successfully imported
    if (pModule != NULL) {
        PyObject *pFunc = PyObject_GetAttrString(pModule, "numpy_function");

        // Check if the function is callable
        if (pFunc && PyCallable_Check(pFunc)) {
            // Create a sample 2D list (matrix) to pass as an argument
            PyObject *matrix = PyList_New(2);

            PyObject *row1 = PyList_New(3);
            PyList_SetItem(row1, 0, PyLong_FromLong(1));
            PyList_SetItem(row1, 1, PyLong_FromLong(2));
            PyList_SetItem(row1, 2, PyLong_FromLong(3));
            PyList_SetItem(matrix, 0, row1);

            PyObject *row2 = PyList_New(3);
            PyList_SetItem(row2, 0, PyLong_FromLong(4));
            PyList_SetItem(row2, 1, PyLong_FromLong(5));
            PyList_SetItem(row2, 2, PyLong_FromLong(6));
            PyList_SetItem(matrix, 1, row2);

            // Create a tuple to hold the arguments
            PyObject *pArgs = PyTuple_Pack(1, matrix);

            // Call the function with arguments
            PyObject *pValue = PyObject_CallObject(pFunc, pArgs);

            // Check for a valid return value from the function
            if (pValue != NULL) {
                printf("Result of call: %ld\n", PyLong_AsLong(pValue));
                Py_DECREF(pValue);
            } else {
                PyErr_Print();
                fprintf(stderr, "Call failed\n");
            }

            // Clean up
            Py_DECREF(pArgs);
            Py_DECREF(matrix);
        } else {
            if (PyErr_Occurred()) PyErr_Print();
            fprintf(stderr, "Cannot find function 'numpy_function'\n");
        }
        Py_XDECREF(pFunc);
        Py_DECREF(pModule);
    } else {
        PyErr_Print();
        fprintf(stderr, "Failed to load \"%s\"\n", "numpy_example");
    }

    // Finalize the Python Interpreter
    Py_Finalize();

    // Confirm end of process
    printf("Process finished successfully.\n");

    return 0;
}