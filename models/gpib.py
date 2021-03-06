# coding: utf8

def sweep_up_btn(field, value):
    return INPUT(_name=field.name,
                 _type='button',
                 _id="%s_%s" % (field._tablename, field.name),
                 _class=field.type,
                 _value='UP',
                 requires=field.requires)
                 
def sweep_down_btn(field, value):
    return INPUT(_name=field.name,
                 _type='button',
                 _id="%s_%s" % (field._tablename, field.name),
                 _class=field.type,
                 _value='DOWN',
                 requires=field.requires)
                 
def calc_sweep_btn(field, value):
    return INPUT(_name=field.name,
                 _type='button',
                 _id="%s_%s" % (field._tablename, field.name),
                 _class=field.type,
                 _value='Calculate',
                 requires=field.requires)                  

db.define_table('gpib',
    Field('device',requires=IS_IN_SET(['Glass','Medium'])),
    Field('run_type',requires=IS_IN_SET(['PSD','Sweep'])),
    Field('multiple_runs','boolean'),
    Field('number_runs','integer'),
    Field('frequency_start',requires=IS_NOT_EMPTY()),
    Field('frequency_range',requires=IS_NOT_EMPTY()),
    Field('source_level',requires=IS_NOT_EMPTY()),
    Field('source_offset',requires=IS_NOT_EMPTY()),
    Field('sweep_rate'),
    Field('calculate_sweep_rate',widget=calc_sweep_btn),
    Field('sweep_up',widget=sweep_up_btn),
    Field('sweep_down',widget=sweep_down_btn))

db.gpib.frequency_range.label = 'Frequency Range/Step Size'
db.gpib.number_runs.label = 'Number of Runs'

db.gpib.source_level.default = '1'
db.gpib.source_offset.default = '1'

db.define_table('results',
    Field('time','datetime'),
    Field('device'),
    Field('run_type'),
    Field('start_frequency'),
    Field('frequency_range'),
    Field('source_level'),
    Field('source_offset'),
	Field('sweep_rate'))

db.define_table('data',
    Field('x'),
    Field('yr'),
    Field('ya'),
    Field('plot_link'))


#db.data.id.represent=lambda id: SPAN(A(IMG(_src=URL(r=request,c='default',f='plot_thumb',args=id),_href=URL(r=request,c='default',f='plot_data',args=id))))
db.data.id.represent=lambda id: SPAN(A(IMG(_src=URL(r=request,c='default',f='plot_thumb',args=id),_width="80px",_height="64px",_title=request.vars.start_frequency),
      _href=URL(r=request,c='default',f='plot_data',args=id),_rel="lightbox[gpib]"))

db.data.id.label='Data Plot'
