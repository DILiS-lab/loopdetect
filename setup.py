import setuptools 

# read the contents of your README file
import pathlib


# The directory containing this file
HERE = pathlib.Path(__file__).parent

# The text of the README file
README = (HERE / "README.md").read_text()


setuptools.setup(
    name='loopdetect',
    version='0.1.0',    
    description='A Python package for feedback loop detection in ODE models',
    url='https://gitlab.com/kabaum/loopdetect', 
    author='Katharina Baum',
    author_email='katharina.baum@hpi.de',
    license='BSD',
    packages=setuptools.find_packages(),
    package_data={'loopdetect': ['data/*.tsv']},
    install_requires=['pandas',  # for output dataframe
                      'numpy',  # for many numerics
                      'numdifftools', # for Jacobian computation
                      #'more-itertools', #we only need itertools which should come with base Python
                      'networkx',       # for path detection
                      'setuptools'       # pkg_resources needs setuptools             
                      ],
    #to adapt
    classifiers=[
     #   'Development Status :: 1 - Planning',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: BSD License',  
        #'Operating System :: OS Independent',        
        #'Programming Language :: Python :: 2',
        #'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.8',
        #'Programming Language :: Python :: 3.5',
    ],
    long_description=README, #insert README text as long description
    long_description_content_type='text/markdown', #parsing as markdown input
    include_package_data=True, #include data for example computations
    project_urls={  # Optional
        'Documentation': 'https://kabaum.gitlab.io/loopdetect',
        #'Source': 'https://kabaum.gitlab.io/loopdetect',
    },
)
