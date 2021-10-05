from setuptools import setup, find_packages

setup(name='tadasv',
      version='0.0.0',
      description='Annotation of TADs and SV',
      long_description='This package allows to annotate SVs with annotations expressing their functional impact. Based on the functional annotation, a classifier can be trained and used to prioritze SVs.',
      classifiers=[
          'Programming Language :: Python :: 3.6',
          ],
      keywords='SV Pathogencity Annotation TADs',
      url='https://github.com/jakob-he/TADA',
      author='Jakob Hertzberg',
      author_email='hertzber@molgen.mpg.de',
      test_suite='tests',
      packages=find_packages(),
      entry_points={
        'console_scripts': [
            'annotate_tads = tadasv.annotate_tads:main',
            'annotate_svs = tadasv.annotate_svs:main',
            'classification_run = tadasv.classification_run:main',
            'predict_pathogenicity = tadasv.predict_pathogenicity:main'
        ]
    },
      install_requires=[
        'scipy>=1.5',
        'numpy>=1.19',
        'pandas>=1.1',
        'matplotlib>=3.3',
        'seaborn>=0.11',
        'pyyaml>=5.4',
        'networkx>=2.5',
        'scikit-learn>=0.24',
        'shap>=0.39.0',
      ],
      include_package_data=True)
