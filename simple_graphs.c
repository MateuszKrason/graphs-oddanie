#define BYTE_TO_BINARY_PATTERN "%c%c%c%c%c%c%c%c"
#define BYTE_TO_BINARY(byte)  \
  (byte & 0x80 ? '1' : '0'), \
  (byte & 0x40 ? '1' : '0'), \
  (byte & 0x20 ? '1' : '0'), \
  (byte & 0x10 ? '1' : '0'), \
  (byte & 0x08 ? '1' : '0'), \
  (byte & 0x04 ? '1' : '0'), \
  (byte & 0x02 ? '1' : '0'), \
  (byte & 0x01 ? '1' : '0')

#include <Python.h>
#include "structmember.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

typedef struct {
    PyObject_HEAD
    short vcs;
    short *edges;
} Adjacency;

static int Adjacency_init(Adjacency *self, PyObject *args, PyObject *kwargs) { //works
    static char *keywords[] = {"text", NULL};
    char *txt = "?";

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "|s", keywords, &txt)) {
        return -1;
    }

    self->edges = malloc(16 * sizeof(short));
    for (int i = 0; i < 16; i++) {
        self->edges[i] = 0x0000;
    }

    int vertices_count = 0;
    for (int i = 0; i < txt[0] - 63; i++) {
        self->vcs = self->vcs << 1;
        self->vcs = self->vcs | 0x0001;
        vertices_count++;
    }

    int var1 = 0;
    int var2 = 1;
    int var3 = 0;

    for (int v = 1; v < vertices_count; v++) {
        for (int u = 0; u < v; u++) {
            if (var1 == 0) {
                var3 = txt[var2] - 63;
                var2++;
                var1 = 6;
            }
            var1--;
            if ((var3 & (1 << var1)) != 0) {
                self->edges[u] = self->edges[u] | (0x0001 << v);
                self->edges[v] = self->edges[v] | (0x0001 << u);
            }
        }
    }

    return 0;
}

static PyObject *Adjacency_impl(PyTypeObject *type, PyObject *args, PyObject *kwds) { //works
    Adjacency *self;
    self = (Adjacency *) type->tp_alloc(type, 0);
    if (self != NULL) {
        self->vcs = 0x0000;
        self->edges = NULL;
    }
    return (PyObject *) self;
}

static void AdjacencyMatrix_unalloc(Adjacency *self) {
    free(self->edges);
    Py_TYPE(self)->tp_free((PyObject *) self);
}

static PyObject *vertices(Adjacency *self) { //corect
    PyObject *vertices_set = PySet_New(NULL);
    if (!vertices_set) {
        return NULL;
    }

    for (int i = 0; i < 16; i++) {
        if (((self->vcs >> i) & 0x0001) == 1) {
            PyObject *item = PyLong_FromLong(i);
            PySet_Add(vertices_set, item);
            Py_DECREF(item);
        }
    }

    return vertices_set;
}

static PyObject *number_of_vertices(Adjacency *self) { //correct
    short count = 0;
    short vcs = self->vcs;
    for (int i = 0; i < 16; i++) {
        if ((vcs & 0x0001) == 1) {
            count++;
        }
        vcs = vcs >> 1;
    }
    long ret = count;
    return PyLong_FromLong(ret);
}

static PyObject *vertex_degree(Adjacency *self, PyObject *args) { //correct
    int v;

    if (args != NULL) {
        PyArg_ParseTuple(args, "i", &v);
    }

    short count = 0;
    short edges = self->edges[v];
    for (int i = 0; i < 16; i++) {
        if ((edges & 0x0001) == 1) {
            count++;
        }
        edges = edges >> 1;
    }

    long ret = (long) count;
    return PyLong_FromLong(ret);
}

static PyObject *vertex_neighbors(Adjacency *self, PyObject *args) { //correct
    int v;

    if (args != NULL) {
        PyArg_ParseTuple(args, "i", &v);
    }
    PyObject *neighbors_set = PySet_New(NULL);

    short edges = self->edges[v];
    for (int i = 0; i < 16; i++) {
        if ((edges & 0x0001) == 1) {
            PyObject *item = PyLong_FromLong(i);
            PySet_Add(neighbors_set, item);
            Py_DECREF(item);
        }
        edges = edges >> 1;
    }

    return neighbors_set;
}

static PyObject *add_vertex(Adjacency *self, PyObject *args) { //correct
    int v;

    if (args != NULL) {
        PyArg_ParseTuple(args, "i", &v);
    }

    short tmp = (0x0001 << v);
    self->vcs = self->vcs | tmp;
    return PyBool_FromLong(1);
}

static PyObject *delete_vertex(Adjacency *self, PyObject *args) { //correct
    int v;

    if (args != NULL) {
        PyArg_ParseTuple(args, "i", &v);
    }

    self->edges[v] = 0x0000;

    short tmp = 0x0001 << v;
    tmp = ~tmp;

    for (int i = 0; i < 16; i++) {
        self->edges[i] = self->edges[i] & tmp;
    }
    self->vcs = self->vcs & tmp;
    return PyBool_FromLong(1);
}

static PyObject *edges(Adjacency *self) { //correct
    PyObject *edges_set = PySet_New(NULL);
    if (!edges_set) {
        return NULL;
    }

    for (int j = 0; j < 16; j++) {
        short edges = self->edges[j];
        for (int i = 0; i < 16; i++) {
            if ((edges & 0x0001) == 1) {
                PyObject *edge = PyTuple_New(2);
                PyTuple_SetItem(edge, 0, PyLong_FromLong(fmin(i, j)));
                PyTuple_SetItem(edge, 1, PyLong_FromLong(fmax(i, j)));
                PySet_Add(edges_set, edge);
                Py_DECREF(edge);
            }
            edges = edges >> 1;
        }
    }

    return edges_set;
}

static PyObject *number_of_edges(Adjacency *self) { //correct
    short count = 0;
    for (int j = 0; j < 16; j++) {
        short edges = self->edges[j];
        for (int i = 0; i < 16; i++) {
            if ((edges & 0x0001) == 1) {
                count++;
            }
            edges = edges >> 1;
        }
    }
    long ret = count / 2;
    return PyLong_FromLong(ret);
}

static PyObject *is_edge(Adjacency *self, PyObject *args) { //correct
    int v, u;

    if (args != NULL) {
        PyArg_ParseTuple(args, "ii", &v, &u);
    }

    short edges_v = self->edges[v];
    edges_v = edges_v >> u;
    edges_v = edges_v & 0x0001;

    return PyBool_FromLong(edges_v);
}

static PyObject *add_edge(Adjacency *self, PyObject *args) { //correct
    int v, u;

    if (args != NULL) {
        PyArg_ParseTuple(args, "ii", &v, &u);
    }

    if (v != u) {
        short update_u = (0x0001 << v);
        short update_v = (0x0001 << u);
        self->edges[v] = self->edges[v] | update_v;
        self->edges[u] = self->edges[u] | update_u;
    }
    return PyBool_FromLong(1);
}

static PyObject *delete_edge(Adjacency *self, PyObject *args) { //correct
    int v, u;

    if (args != NULL) {
        PyArg_ParseTuple(args, "ii", &v, &u);
    }

    short update_u = ~(0x0001 << v);
    short update_v = ~(0x0001 << u);
    self->edges[v] = self->edges[v] & update_v;
    self->edges[u] = self->edges[u] & update_u;
    return PyBool_FromLong(1);
}
static int getN(Adjacency *self) {
    short n = 0;
    short vcs = self->vcs;
    for (int i = 0; i < 16; i++) {
        if ((vcs & 0x0001) == 1) {
            n++;
        }
        vcs = vcs >> 1;
    }
    return n;
}

static void mark_component(Adjacency* self, int vertex, int* cs, int mark)
{
    cs[vertex] = mark;
    int n = getN(self);
    for (int v = 0; v < n; v++) {
        if (cs[v] == 0 && PyObject_IsTrue(is_edge(self, Py_BuildValue("(ii)", vertex, v)))) {
            mark_component(self, v, cs, mark);
        }
    }
}

static PyObject *connected_components(Adjacency *self){ //connected
    if (self->vcs == 0x00) {
        return PyLong_FromLong(0);
    }
    int n = getN(self);
    int cs[n];
    memset(cs, 0, n*sizeof(int));
    int c = 0;
    for (int v = 0; v < n; v++) {
        if (cs[v] == 0) {
            c += 1;
            mark_component(self, v, cs, c);
        }
    }
    return PyLong_FromLong(c);
}


static PyMemberDef Adjacency_members[] = {
        {"vertices", T_SHORT, offsetof(Adjacency, vertices), 0, PyDoc_STR("vertices of the graph")},
        {"edges",    T_SHORT, offsetof(Adjacency, edges),    0, PyDoc_STR("edges of the graph")},
        {NULL}
};

static PyMethodDef Adjacency_methods[] = {
        {"number_of_vertices", (PyCFunction) number_of_vertices, METH_NOARGS},
        {"vertices",           (PyCFunction) vertices,           METH_NOARGS},
        {"vertex_degree",      (PyCFunction) vertex_degree,      METH_VARARGS},
        {"vertex_neighbors",   (PyCFunction) vertex_neighbors,   METH_VARARGS},
        {"add_vertex",         (PyCFunction) add_vertex,         METH_VARARGS},
        {"delete_vertex",      (PyCFunction) delete_vertex,      METH_VARARGS},
        {"number_of_edges",    (PyCFunction) number_of_edges,    METH_NOARGS},
        {"is_edge",            (PyCFunction) is_edge,            METH_VARARGS},
        {"add_edge",           (PyCFunction) add_edge,           METH_VARARGS},
        {"delete_edge",        (PyCFunction) delete_edge,        METH_VARARGS},
        {"edges",              (PyCFunction) edges,              METH_NOARGS},
        {"connected_components",       (PyCFunction) connected_components,       METH_NOARGS},
        {NULL,                 NULL}
};

static PyTypeObject AdjacencyMatrixType = {
        PyVarObject_HEAD_INIT(NULL, 0)
        "_simple_graphs.AdjacencyMatrix",
        sizeof(Adjacency),
        0,
        (destructor)AdjacencyMatrix_unalloc,
        0,                                  
        0,                                  
        0,                                  
        0,                                  
        0,                                  
        0,                                  
        0,                                  
        0,                                  
        0,                                  
        0,                                  
        0,                                  
        0,                                  
        0,                                  
        0,                                  
        Py_TPFLAGS_DEFAULT,                 
        0,                                  
        0,                                  
        0,                                  
        0,                                  
        0,                                  
        0,                                  
        0,                                  
        Adjacency_methods,                  
        Adjacency_members,            
        0,                                  
        0,                                  
        0,                                  
        0,                                  
        0,                                  
        0,                                  
        (initproc)Adjacency_init,     
        0,                                  
        Adjacency_impl,               
};

static struct PyModuleDef graphmodule = {
        PyModuleDef_HEAD_INIT,
        "_simple_graphs",
        NULL,
        -1
};

PyMODINIT_FUNC PyInit_simple_graphs(void) {
    PyObject *m;
    if (PyType_Ready(&AdjacencyMatrixType) < 0)
        return NULL;

    m = PyModule_Create(&graphmodule);
    if (m == NULL)
        return NULL;

    Py_INCREF(&AdjacencyMatrixType);
    if (PyModule_AddObject(m, "AdjacencyMatrix", (PyObject * ) & AdjacencyMatrixType) < 0) {
        Py_DECREF(&AdjacencyMatrixType);
        Py_DECREF(m);
        return NULL;
    }

    return m;
}