from distutils.core import setup
from Cython.Build import cythonize

setup(
    name        = "EnergyModel",
    version     = '0.0.0',
    author      = 'Ruben Tikidji-Hamburyan',
    author_email= 'rath@gwu.edu',
    py_modules  = [".EnergyModel"],
    ext_modules = cythonize("EnergyChaser.pyx"),
    description = """
/***********************************************************************************************************\

This setup script builds EnergyChaser library associated with paper:                                                        
 Ruben Tikidji-Hamburyan , Tarek El-Ghazawi , Jason Triplett
    Novel Models of Visual Topographic Map Alignment into Superior Colliculus

 Copyright: Ruben Tikidji-Hamburyan <rath@gwu.edu> Apr.2016 - Sep.2016

\************************************************************************************************************/    
    """
    
)

