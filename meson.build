project('ephe', 'cpp', license: 'GPL3+', version : '2.0')

intldep = dependency('intl', required: false)

conf_data = configuration_data()
conf_data.set('HAVE_LIBINTL', intldep.found())
configure_file(output : 'config.h', configuration: conf_data)
confinc = include_directories('.')

executable('ephe',
  ['main.cpp', 'vect.cpp', 'body.cpp', 'orbit.cpp',
  'observ.cpp', 'instant.cpp', 'ephe.cpp'],
  dependencies: intldep,
  include_directories: confinc,
  install: true)

