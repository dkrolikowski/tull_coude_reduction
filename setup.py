from setuptools import setup

setup( name = 'tull_coude_reduction',
       description = 'Tull coude spectrograph reduction and analysis pipeline',
       author = 'Danny Krolikowski',
       author_email = 'krolikowski@arizon.edu',
       url = 'https://github.com/krolikowski/tull_coude_reduction',
       license = 'MIT',
       classifiers = [ 'Development Status :: 3 - Alpha',
                       'Intended Audience :: Science/Research',
                       'License :: OSI Approved :: MIT License',
                       'Programming Language :: Python :: 3',
                       'Topic :: Scientific/Engineering :: Astronomy',
                       ],
       project_urls = { 'Documentation': 'https://tull-coude-reduction.readthedocs.io/en/latest/'
                       },
       packages = [ 'tull_coude_reduction', 'tull_coude_reduction.modules' ],
       python_requires = '>=3',
       install_requires = [ 'numpy',
                            'matplotlib',
                            'scipy',
                            'pandas',
                            'astropy==5.2',
                            'astroscrappy',
                            'saphires',
                            'barycorrpy',
                            'tqdm',
                            'yaml'
                            ]
       )