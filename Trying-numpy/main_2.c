//
// Created by Mofe Fagade on 9/19/24.
//

// main_2.c
/*START_USER_CODE*/
#include <Python.h>
//#include <numpy/arrayobject.h>
#include <stdio.h>

    int main(int argc, char *argv[]) {
        Py_Initialize();

        if (!Py_IsInitialized()) {
            fprintf(stderr, "Python initialization failed\n");
            return 1;
        }

        double rbInitial = 0.0215e-6;
        double rbDotInitial = 0.1; //adjustment
        double CL_ = 1300.0;
        double Pv_ = 1160000;
        double pAmbient = 116000;
        double sigma = 0.25;
        double n_ = 1.0;
        double mu = 0.24e-3;
        double kappa_ = 8.0e-6;
        double time = 1.5e-7; //adjustment
        //double rb = 0.0;
        double rho = 595.59;

        PyObject* pModule = NULL;
        PyObject* pFunc = NULL;

        printf("this is an update\n");
        PyRun_SimpleString("import sys");
        printf("this is an updated sys\n");
        PyRun_SimpleString("sys.path.append(\".\")");
        printf("this is an updated sys.path\n");

        printf("Importing KMEquation module\n");
        pModule = PyImport_ImportModule("KMEquation");
        if (!pModule) {
            fprintf(stderr, "Module or script file not found\n");
            if (PyErr_Occurred()) PyErr_Print();
            Py_Finalize();
            return 1;
        }

        printf("Finding the solve function in python\n");
        pFunc = PyObject_GetAttrString(pModule, "solve");
        if (!pFunc || !PyCallable_Check(pFunc)) {
            fprintf(stderr, "No callable function 'solve' found\n");
            if (PyErr_Occurred()) PyErr_Print();
            Py_Finalize();
            return 1;
        }

        double initial_conditions[] = {rbInitial, rbDotInitial};
        PyObject* pyInitialConditions = Py_BuildValue("(dd)", initial_conditions[0], initial_conditions[1]);
        PyObject* pyList = PyList_New(2);
        PyList_SetItem(pyList, 0, PyFloat_FromDouble(0.0));
        PyList_SetItem(pyList, 1, PyFloat_FromDouble(time));

        if (!pyInitialConditions || !pyList) {
            fprintf(stderr, "Argument preparation failed\n");
            if (PyErr_Occurred()) PyErr_Print();
            Py_Finalize();
            return 1;
        }

        PyObject* pyArgs = PyTuple_Pack(13, pyInitialConditions, PyFloat_FromDouble(rbInitial), PyFloat_FromDouble(CL_), PyFloat_FromDouble(Pv_), PyFloat_FromDouble(pAmbient), PyFloat_FromDouble(pAmbient), PyFloat_FromDouble(sigma), PyFloat_FromDouble(n_), PyFloat_FromDouble(mu), PyFloat_FromDouble(kappa_), PyFloat_FromDouble(rho), pyList, PyFloat_FromDouble(time));

        if (!pyArgs) {
            fprintf(stderr, "Packing arguments failed\n");
            if (PyErr_Occurred()) PyErr_Print();
            Py_Finalize();
            return 1;
        }

        printf("Calling solve function\n");
        PyObject* pReturn;
        pReturn = PyObject_CallObject(pFunc, pyArgs);

        if (!pReturn) {
            fprintf(stderr, "Function call failed\n");
            if (PyErr_Occurred()) PyErr_Print();
        } else {
            printf("Returned value type: %s\n", Py_TYPE(pReturn)->tp_name);

            if (PyList_Check(pReturn)) {
                Py_ssize_t size = PyList_Size(pReturn);
                for (Py_ssize_t i = 0; i < size; ++i) {
                    PyObject* item = PyList_GetItem(pReturn, i);
                    // Assume all items are doubles
                    if (PyFloat_Check(item)) {
                        double value = PyFloat_AsDouble(item);
                        printf("R_value[%ld]: %f\n", i, value);
                    } else {
                        fprintf(stderr, "Non-float item in R_values list\n");
                    }
                }
            } else {
                fprintf(stderr, "Returned value is not a list\n");
            }

            double rb;
            if (PyArg_Parse(pReturn, "d", &rb)) {
                printf("rb: %f\n", rb);
            }
        }

        Py_XDECREF(pyArgs);
        Py_XDECREF(pyInitialConditions);
        Py_XDECREF(pyList);
        Py_XDECREF(pReturn);
        Py_XDECREF(pModule);
        Py_XDECREF(pFunc);

        Py_Finalize();

        return 0;
    }

