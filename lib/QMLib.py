import os
import numpy as np

import ctypes

class QMLib(object):
    def __init__(self):
        self.qm = np.ctypeslib.load_library("./libfort.so", os.path.dirname(os.path.abspath(__file__)))

    def set_time(self, time):
        self.qm.set_time.argtypes = [ctypes.POINTER(ctypes.c_double)]
        self.qm.set_time.restype = ctypes.c_void_p
        time = ctypes.byref(ctypes.c_double(time))
        self.qm.set_time(time)

    def get_time(self):
        self.qm.get_time.argtypes = [ctypes.POINTER(ctypes.c_double)]
        self.qm.get_time.restype = ctypes.c_void_p
        time = ctypes.c_double(0.0)
        self.qm.get_time(ctypes.byref(time))
        return(time.value)
        
    def set_dx(self, dx):
        self.qm.set_dx.argtypes = [ctypes.POINTER(ctypes.c_double)]
        self.qm.set_dx.restype = ctypes.c_void_p
        dx = ctypes.byref(ctypes.c_double(dx))
        self.qm.set_dx(dx)

    def get_dx(self):
        self.qm.get_dx.argtypes = [ctypes.POINTER(ctypes.c_double)]
        self.qm.get_dx.restype = ctypes.c_void_p
        dx = ctypes.c_double(0.0)
        self.qm.get_dx(ctypes.byref(dx))
        return(dx.value)

    def set_dy(self, dy):
        self.qm.set_dy.argtypes = [ctypes.POINTER(ctypes.c_double)]
        self.qm.set_dy.restype = ctypes.c_void_p
        dy = ctypes.byref(ctypes.c_double(dy))
        self.qm.set_dy(dy)

    def get_dy(self):
        self.qm.get_dy.argtypes = [ctypes.POINTER(ctypes.c_double)]
        self.qm.get_dy.restype = ctypes.c_void_p
        dy = ctypes.c_double(0.0)
        self.qm.get_dy(ctypes.byref(dy))
        return(dy.value)
    
    def set_nx(self, nx):
        self.qm.set_nx.argtypes = [ctypes.POINTER(ctypes.c_int32)]
        self.qm.set_nx.restype = ctypes.c_void_p
        nx = ctypes.byref(ctypes.c_int32(nx))
        self.qm.set_nx(nx)

    def get_nx(self):
        self.qm.get_nx.argtypes = [ctypes.POINTER(ctypes.c_int32)]
        self.qm.get_nx.restype = ctypes.c_void_p
        nx = ctypes.c_int32(0)
        self.qm.get_nx(ctypes.byref(nx))
        return(nx.value)

    def set_ny(self, ny):
        self.qm.set_ny.argtypes = [ctypes.POINTER(ctypes.c_int32)]
        self.qm.set_ny.restype = ctypes.c_void_p
        ny = ctypes.byref(ctypes.c_int32(ny))
        self.qm.set_ny(ny)

    def get_ny(self):
        self.qm.get_ny.argtypes = [ctypes.POINTER(ctypes.c_int32)]
        self.qm.get_ny.restype = ctypes.c_void_p
        ny = ctypes.c_int32(0)
        self.qm.get_ny(ctypes.byref(ny))
        return(ny.value)

    def set_psi2d(self, psi2d_real, psi2d_imag):
        self.qm.set_psi2d.argtypes = [np.ctypeslib.ndpointer(dtype=np.float64),
                                      np.ctypeslib.ndpointer(dtype=np.float64)]
        self.qm.set_psi2d.restype = ctypes.c_void_p
        self.qm.set_psi2d(psi2d_real, psi2d_imag)
        
    def get_psi2d(self, psi2d_real, psi2d_imag):
        self.qm.get_psi2d.argtypes = [np.ctypeslib.ndpointer(dtype=np.float64),
                                      np.ctypeslib.ndpointer(dtype=np.float64)]
        self.qm.get_psi2d.restype = ctypes.c_void_p
        self.qm.get_psi2d(psi2d_real, psi2d_imag)

    def set_v2d(self, v2d_real, v2d_imag):
        self.qm.set_v2d.argtypes = [np.ctypeslib.ndpointer(dtype=np.float64),
                                      np.ctypeslib.ndpointer(dtype=np.float64)]
        self.qm.set_v2d.restype = ctypes.c_void_p
        self.qm.set_v2d(v2d_real, v2d_imag)
        
    def get_v2d(self, v2d_real, v2d_imag):
        self.qm.get_v2d.argtypes = [np.ctypeslib.ndpointer(dtype=np.float64),
                                      np.ctypeslib.ndpointer(dtype=np.float64)]
        self.qm.get_v2d.restype = ctypes.c_void_p
        self.qm.get_v2d(v2d_real, v2d_imag)
        
    def set_field_x(self, field_x):
        self.qm.set_field_x.argtypes = [np.ctypeslib.ndpointer(dtype=np.float64)]
        self.qm.set_field_x.restype = ctypes.c_void_p
        self.qm.set_field_x(field_x)
        
    def get_field_x(self, field_x):
        self.qm.get_field_x.argtypes = [np.ctypeslib.ndpointer(dtype=np.float64)]
        self.qm.get_field_x.restype = ctypes.c_void_p
        self.qm.get_field_x(field_x)

    def set_field_y(self, field_y):
        self.qm.set_field_y.argtypes = [np.ctypeslib.ndpointer(dtype=np.float64)]
        self.qm.set_field_y.restype = ctypes.c_void_p
        self.qm.set_field_y(field_y)
        
    def get_field_y(self, field_y):
        self.qm.get_field_y.argtypes = [np.ctypeslib.ndpointer(dtype=np.float64)]
        self.qm.get_field_y.restype = ctypes.c_void_p
        self.qm.get_field_y(field_y)

    def init_2d(self):
        self.qm.init_2d.argtypes = []
        self.qm.init_2d.restype = ctypes.c_void_p
        self.qm.init_2d()
        
    def put_gaussian(self, x, y, kx, ky, c):
        self.qm.put_gaussian.argtypes = [ctypes.POINTER(ctypes.c_double),
                                         ctypes.POINTER(ctypes.c_double),
                                         ctypes.POINTER(ctypes.c_double),
                                         ctypes.POINTER(ctypes.c_double),
                                         ctypes.POINTER(ctypes.c_double)]
        self.qm.put_gaussian.restype = ctypes.c_void_p
        x = ctypes.byref(ctypes.c_double(x))
        y = ctypes.byref(ctypes.c_double(y))
        kx = ctypes.byref(ctypes.c_double(kx))
        ky = ctypes.byref(ctypes.c_double(ky))
        c = ctypes.byref(ctypes.c_double(c))
        self.qm.put_gaussian(x, y, kx, ky, c)

    def set_v2d_box(self, sx, sy, ex, ey, Re_c, Im_c):
        self.qm.set_v2d_box.argtypes = [ctypes.POINTER(ctypes.c_double),
                                        ctypes.POINTER(ctypes.c_double),
                                        ctypes.POINTER(ctypes.c_double),
                                        ctypes.POINTER(ctypes.c_double),
                                        ctypes.POINTER(ctypes.c_double),
                                        ctypes.POINTER(ctypes.c_double)]
        self.qm.set_v2d_box.restype = ctypes.c_void_p
        sx = ctypes.byref(ctypes.c_double(sx))
        sy = ctypes.byref(ctypes.c_double(sy))
        ex = ctypes.byref(ctypes.c_double(ex))
        ey = ctypes.byref(ctypes.c_double(ey))
        Re_c = ctypes.byref(ctypes.c_double(Re_c))
        Im_c = ctypes.byref(ctypes.c_double(Im_c))
        self.qm.set_v2d_box(sx, sy, ex, ey, Re_c, Im_c)

    def step_2d(self, n):
        self.qm.step_2d.argtypes = [ctypes.POINTER(ctypes.c_int32)]
        self.qm.step_2d.restype = ctypes.c_void_p
        n = ctypes.byref(ctypes.c_int32(n))
        self.qm.step_2d(n)

