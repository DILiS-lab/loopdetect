from setuptools import setup

# read the contents of your README file
from os import path
this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()


setup(
    name='loopdetect',
    version='0.1.0',    
    description='A Python package for feedback loop detection in ODE models',
    url='https://gitlab.com/kabaum/loopdetect', 
    author='Katharina Baum',
    author_email='katharina.baum@hpi.de',
    license='GPL-3',
    packages=['loopdetect'],
    install_requires=['pandas', 
                      'numpy', 
                      'numdifftools',
                      'itertools',
                      'networkx',
                      'pkg_references'                    
                      ],
    #to adapt
    classifiers=[
     #   'Development Status :: 1 - Planning',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: BSD License',  
        #'Operating System :: POSIX :: Linux',        
        #'Programming Language :: Python :: 2',
        #'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.8',
        #'Programming Language :: Python :: 3.5',
    ],
    long_description=long_description, #insert README text as long description
    long_description_content_type='text/markdown', #parsing as markdown input
    include_package_data=True #include data for example computations
)
