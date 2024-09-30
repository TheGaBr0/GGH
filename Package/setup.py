from distutils.core import setup
setup(
  name = 'GGH_crypto',
  packages = ['GGH_crypto'],
  version = '1',
  license='MIT',
  description = 'GGH_crypto is a Python package for lattice-based cryptography, focusing on GGH and GGH-HNF implementations.',
  author = 'Gabriele Bottani',
  author_email = 'gbotani19@gmail.com',
  url = 'https://github.com/TheGaBr0/GGH',
  download_url = 'https://github.com/TheGaBr0/GGH/archive/refs/tags/v1.0.tar.gz',
  keywords = ['GGH', 'GGH-HNF', 'GGH_CRYPTO', 'Lattice', 'LLL', 'BKZ', 'Lattice-based-cryptography'],
  install_requires=[
      'matplotlib==3.9.2',
      'numpy==2.1.1',
      'python_flint==0.6.0',
      'sympy==1.12.1',
  ],
  classifiers=[
    'Development Status :: 5 - Production/Stable',
    'Intended Audience :: Developers',
    'Topic :: Software Development :: Build Tools',
    'License :: OSI Approved :: MIT License',
    'Programming Language :: Python :: 3.10',
    'Programming Language :: Python :: 3.11',
    'Programming Language :: Python :: 3.12'
  ],
)