#include <Python.h>
#include <math.h>

static PyObject * errorObject; 

#define SET_ERROR (message) \
   {                       \
       PyErr_SetString(errorObject, message); \
   }                   
#define ON_ERROR (message) \
   {                       \
       PyErr_SetString(errorObject, message); \
       return NULL;                          \
   }


static
PyObject * GetTrans (PyObject * self, PyObject * args) {
          int Lx, Ly, Nt; 
          int wrap; 
	  int i;
          const int maxnt = 10000; 
          double energy, flux;
          double *transvals;
          double check; 
          char *current_dir; 
          char *gauge_dir; 
          double *V; 
          PyObject *retval; 
          PyObject *tval_list;
          PyObject *vlist; 
          int iV, jV; 
   
          double gettrans_ (char *current_dir, char *gauge_dir, 
			    double *t, int *nt, int *Lx, int *Ly,
			    double *V, 
			    double *E, double *F, int *w); 
   
          if (!PyArg_ParseTuple (args, "ssiiddiO", 
				 &current_dir, &gauge_dir, 
				 &Lx, &Ly, &energy, &flux, &wrap, 
				 &vlist)) {
	       PyErr_SetString(errorObject, "invalid arguments to GetTrans");
	       return NULL; 
	  }
          
          if (!PySequence_Check(vlist)) {
	       PyErr_SetString(errorObject, "V is not a sequence"); 
	       return NULL; 
	  }
   
          V = PyMem_New (double, Lx * Ly);
          for (iV = 0; iV < Lx; iV++) {
	      PyObject * row = PySequence_GetItem(vlist, iV); 
	      if (!PySequence_Check(row)) {
	           PyErr_SetString(errorObject, "V[i] not a sequence");
		   PyMem_Free(V); 
		   return NULL; 
	      }
	      for (jV = 0; jV < Ly; jV++) {
	          //V[iV * Ly + jV] = 0.0;
		 PyObject * value = PySequence_GetItem(row, jV);
		 double vij; 
		 if (!PyArg_Parse(value, "d", &vij)) {
		     PyErr_SetString(errorObject, "cannot parse V[i][j]"); 
		     PyMem_Free(V); 
		     return NULL; 
		 }
		 V[iV * Ly + jV] = vij; 
	      }
	  }
   
          transvals = PyMem_New (double, maxnt); 
          check = gettrans_ (current_dir, gauge_dir, transvals, &Nt, 
			     &Lx, &Ly, V, &energy, &flux, &wrap); 
          PyMem_Free(V); 
          tval_list = PyList_New (Nt); 
          for (i = 0; i < Nt; i++) {
	       PyList_SetItem (tval_list, i, 
			       Py_BuildValue ("d", transvals[i])); 
	  }
   
          PyMem_Free (transvals); 
   
          retval = Py_BuildValue ("Nd", tval_list, check); 
          return retval; 
}


static struct PyMethodDef module_methods[] = {
   { "gettrans", GetTrans, METH_VARARGS, "Calculates transmission values" }, 
   { NULL,       NULL,     0,            NULL}
};


void init_tmatrix () {
     PyObject *module, *moduleDict; 
     PyMethodDef *def; 
     module = Py_InitModule ("_tmatrix", module_methods); 
     moduleDict = PyModule_GetDict (module); 
     errorObject = Py_BuildValue("s", "_tmatrix.error"); 
     if (PyErr_Occurred()) {
         Py_FatalError ("cannot initialize the module _tmatrix"); 
     }
}
