from setuptools import setup, find_packages


setup(name='TADA',
      version='0.1',
      description='Annotation of TADs and CNVs',
      long_description='This package allows to annotate CNVs with annotations expressing their functional impact. Based on the functional annotation, a classifier can be trained and used to prioritze CNVs.',
      classifiers=[
          'Development Status :: 3 - Alpha',
          'Programming Language :: Python :: 3.7',
          'Topic :: ',
      ],
      keywords='CNV Pathogencity Annotation TADs',
      url='https://github.com/jakob-he/TADA',
      author='Jakob Hertzberg',
      author_email='jakob.hertzberg@gmail.com',
      test_suite='tests',
      packages=find_packages(),
      entry_points={
        'console_scripts': [
            'annotate_tads = annotate_tads:main',
            'annotate_cnvs = annotate_cnvs:main',
            'classification_run = classification_run:main',
            'predict_variants = predict_variants:main'
        ]
    },
      install_requires=[
        'scipy',
        'numpy',
        'pandas',
        'sklearn',
        'imblearn',
        'matplotlib',
        'seaborn',
        'networkx'
      ],
      include_package_data=True,
      zip_safe=False)
