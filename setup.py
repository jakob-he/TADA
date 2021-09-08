from setuptools import setup, find_packages

setup(name='tadacnv',
      version='1.0.1',
      description='Annotation of TADs and CNVs',
      long_description='This package allows to annotate CNVs with annotations expressing their functional impact. Based on the functional annotation, a classifier can be trained and used to prioritze CNVs.',
      classifiers=[
          'Programming Language :: Python :: 3.6',
          ],
      keywords='CNV Pathogencity Annotation TADs',
      url='https://github.com/jakob-he/TADA',
      author='Jakob Hertzberg',
      author_email='jakob.hertzberg@gmail.com',
      test_suite='tests',
      packages=find_packages(),
      entry_points={
        'console_scripts': [
            'annotate_tads = tadacnv.annotate_tads:main',
            'annotate_cnvs = tadacnv.annotate_cnvs:main',
            'classification_run = tadacnv.classification_run:main',
            'predict_variants = tadacnv.predict_variants:main'
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
        'scikit-learn>=0.24'
      ],
      include_package_data=True,
      zip_safe=False)
