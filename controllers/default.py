# -*- coding: utf-8 -*-

#########################################################################
## This is a samples controller
## - index is the default action of any application
## - user is required for authentication and authorization
## - download is for downloading files uploaded in the db (does streaming)
## - call exposes all registered services (none by default)
#########################################################################
import time
import socket
import numpy as np
import scipy as sp
import math
import binascii
import datetime as dt
import pickle
from time import localtime, strftime
from pylab import *

@auth.requires_login()
def index():
    form = SQLFORM(db.gpib,fields=['device','run_type', 'multiple_runs','number_runs','frequency_start','frequency_range','source_level','source_offset','sweep_rate','calculate_sweep_rate','sweep_up','sweep_down'],col3={'frequency_start':'Hz','frequency_range':'Hz','source_level':'V','source_offset':'V','sweep_rate':'Hz/s'})
    if form.accepts(request.vars, session):
        device = request.vars.device

        if device == "Medium":
            HOST = "172.21.7.68"
        else:
            HOST = "172.21.7.69"

            PORT = 1234
            NL = "\r\n"

            SF_pre = "SF "
            FRS_pre = "FRS "
            SRC_pre = "SRLV "
            SROF_pre = "DCOF "
            V_suf = " V;" + NL
            HZ_suf = " HZ;" + NL

            sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM, socket.IPPROTO_TCP)
            sock.settimeout(0.5)
            sock.connect((HOST,PORT))

            sock.send("++mode 1" + NL)
            sock.send("++addr 20" + NL)
            sock.send("++auto 0" + NL)
            sock.send("++eoi 1" + NL)

            run_type = request.vars.run_type

            if run_type == "PSD":

                SF_in = request.vars.frequency_start
                FRS_in = request.vars.frequency_range

                SF_var = binascii.a2b_qp(SF_pre+str(SF_in)+HZ_suf)
                FRS_var = binascii.a2b_qp(FRS_pre+str(FRS_in)+HZ_suf)

                sock.send("SNGL;" + NL)
                time.sleep(0.5)
                sock.send("MSMD;" + NL)
                time.sleep(0.5)
                sock.send("LNRS;" + NL)
                time.sleep(0.5)
                sock.send("MDSP;" + NL)
                time.sleep(0.5)
                sock.send("PSP2;" + NL)
                time.sleep(0.5)
                sock.send("SMES;" + NL)
                time.sleep(0.5)
                sock.send("PSPC;" + NL)
                time.sleep(0.5)
                sock.send("CH2;" + NL)
                time.sleep(0.5)
                sock.send("FREQ;" + NL)
                time.sleep(0.5)
                sock.send(SF_var)
                time.sleep(0.5)
                sock.send(FRS_var)
                time.sleep(0.5)
                sock.send("CORD;" + NL)
                time.sleep(0.5)
                sock.send("MGDB;" + NL)
                time.sleep(0.5)
                sock.send("SCAL;" + NL)
                time.sleep(0.5)
                sock.send("A;" + NL)
                time.sleep(0.5)
                sock.send("YASC;" + NL)
                time.sleep(0.5)
                sock.send("STRT;" + NL)
                time.sleep(1)
                sock.send("AS?" + NL)
                time.sleep(2)
                sock.send("++read eoi" + NL)
                filter_check = sock.recv(64)
                filter_check = float(filter_check)
                time.sleep(10)
                while filter_check == 2.0:
                    time.sleep(10)
                    sock.send("AS?" + NL)
                    time.sleep(3)
                    sock.send("++read eoi" + NL)
                    time.sleep(4)
                    filter_check = sock.recv(64)
                    time.sleep(20)
                    sock.send("PAUS;" + NL)
                    time.sleep(5)

                # Active Data Trace Dump in ASCII
                sock.send("DDAS;" + NL)
                time.sleep(10)
                sock.send("++read eoi" + NL)
                time.sleep(10)
                dat_headerDDAS = sock.recv(32768)
                dat_headerDDAS = [dat_headerDDAS]
                dat_headerDDAS = dat_headerDDAS[0].split('\r\n')
                headerDDAS = dat_headerDDAS[:67]
                dat =  dat_headerDDAS[67:]
                dat.pop()
                for index, item in enumerate(dat):
                    dat[index] = float(item)

                # Instrument State Dump in ASCII
                sock.send("DSAS;" + NL)
                time.sleep(10)
                sock.send("++read eoi" + NL)
                time.sleep(10)
                state_header = sock.recv(32768)
                state_header = [state_header]
                state_header = state_header[0].split('\r\n')
                state_header.pop()

                # Retrieve useful information from the headers before generating plot
                frq_span = float(state_header[80])
                frq_strt = float(state_header[92])
                frq_res = float(state_header[82])
                x_cord = float(headerDDAS[11])

                # Scale the data as need be (currently only scales for dB)
                for i,n in enumerate(dat):
                    if dat[i] == 0:
                        n = 0;
                    else:
                        dat[i] = 10*math.log10(n)

                # Generate Frequency values based on start frequncy, frequncy span and frequency resolution
                x_val = np.arange(frq_strt,frq_strt+frq_span,frq_res)
                x_val = np.append(x_val,x_val[-1]+frq_res)

                # Zero check and array resize
                try:
                    while (dat.index(0)):
                        zero_check = dat.index(0)
                        dat.remove(0)
                        x_val = np.delete(x_val,zero_check)
                except ValueError:
                    pass

            else:

                mruns = request.vars.multiple_runs

                if (mruns == False):

                    sock.send("STRT;" + NL)
                    time.sleep(20)
                    sock.send("SMSD;" + NL)
                    time.sleep(10)
                    sock.send("++read eoi" + NL)
                    time.sleep(2)
                    ms_check = sock.recv(64)

                    while ms_check == '0\r\n':
                        time.sleep(10)
                        sock.send("SMSD;" + NL)
                        time.sleep(10)
                        sock.send("++read eoi" + NL)
                        time.sleep(4)
                        ms_check = sock.recv(64)
                        time.sleep(1)

                    # Active Data Trace Dump in ASCII
                    sock.send("DDAS;" + NL)
                    time.sleep(15)
                    sock.send("++read eoi" + NL)
                    time.sleep(15)
                    sin_dat_header = sock.recv(131072)
                    sin_dat_header = [sin_dat_header]
                    sin_dat_header = sin_dat_header[0].split('\r\n')
                    sin_dat_header.pop()
                    sin_header = sin_dat_header[:67]
                    sin_dat =  sin_dat_header[67:]
                    for index, item in enumerate(sin_dat):
                        sin_dat[index] = float(item)

                    # Create data arrays using the dumped real and imaginary pairs
                    sin_dat = np.array(sin_dat)
                    sin_dat.shape = 801,2
                    sin_dat_real = sin_dat[:,0]
                    sin_dat_imag = sin_dat[:,1]

                    sin_dat_comp = []
                    for i in range(len(sin_dat_real)):
                        sin_dat_comp.append(complex(sin_dat_real[i],sin_dat_imag[i]))

                    sin_dat_abs = []
                    for i in range(len(sin_dat_comp)):
                        sin_dat_abs.append(abs(sin_dat_comp[i]))


                    # Instrument State Dump in ASCII
                    sock.send("DSAS;" + NL)
                    time.sleep(15)
                    sock.send("++read eoi" + NL)
                    time.sleep(15)
                    sin_state = sock.recv(131072)
                    sin_state = [sin_state]
                    sin_state = sin_state[0].split('\r\n')
                    sin_state.pop()

                    # Retrieve useful information from the headers before generating plot
                    frq_strt = float(sin_state[94])
                    frq_end = float(sin_state[95])
                    frq_span = frq_end - frq_strt
                    frq_res = frq_span/801
                    x_cord = float(sin_header[11])
                    swp_rate = float(sin_state[68])
                    swp_time = float(sin_state[67])
                    swp_res = float(sin_state[69])
                    swp_intm = float(sin_state[70])

                    # Generate Frequency values based on start frequncy, frequncy span and frequency resolution
                    x_val_ss = np.arange(frq_strt,frq_strt+frq_span,frq_res)

                    session.run_type = request.vars.run_type

                    if run_type == "PSD":
                        x_psd_pck = pickle.dumps(x_val)
                        x_psd_dat = pickle.dumps(dat)
                    else:
                        x_pck = pickle.dumps(x_val_ss)
                        yr_pck = pickle.dumps(sin_dat_real)
                        ya_pck = pickle.dumps(sin_dat_abs)
                        db.data.insert(x = x_pck,
                                       yr = yr_pck,
                                       ya = ya_pck)

                else:

                    mSF_in = request.vars.frequency_start
                    mStep_size = request.vars.frequency_range
                    mNumber = request.vars.number_runs

                    mSF_array = []
                    mSF_db_array = []

                    while mNumber > 0:
                        mSF_array.append(binascii.a2b_qp(SF_pre+str(mSF_in)+HZ_suf))
                        mSF_db_array.append(int(float(mSF_in)))
                        mSF_in = int(float(mSF_in)) + int(float(mStep_size))
                        mNumber = int(float(mNumber)) - 1

                    for n,item in enumerate(mSF_array):
                        sock.send("FREQ;" + NL)
                        time.sleep(0.5)
                        sock.send(mSF_array[n])
                        time.sleep(0.5)
                        sock.send("STRT;" + NL)
                        time.sleep(20)
                        sock.send("SMSD;" + NL)
                        time.sleep(10)
                        sock.send("++read eoi" + NL)
                        time.sleep(2)
                        ms_check = sock.recv(64)
                    
                        while ms_check == '0\r\n':
                            time.sleep(10)
                            sock.send("SMSD;" + NL)
                            time.sleep(10)
                            sock.send("++read eoi" + NL)
                            time.sleep(4)
                            ms_check = sock.recv(64)
                            time.sleep(1)

                        # Active Data Trace Dump in ASCII
                        sock.send("DDAS;" + NL)
                        time.sleep(15)
                        sock.send("++read eoi" + NL)
                        time.sleep(15)
                        sin_dat_header = sock.recv(131072)
                        sin_dat_header = [sin_dat_header]
                        sin_dat_header = sin_dat_header[0].split('\r\n')
                        sin_dat_header.pop()
                        sin_header = sin_dat_header[:67]
                        sin_dat =  sin_dat_header[67:]
                        for index, item in enumerate(sin_dat):
                            sin_dat[index] = float(item)

                        # Create data arrays using the dumped real and imaginary pairs
                        sin_dat = np.array(sin_dat)
                        sin_dat.shape = 801,2
                        sin_dat_real = sin_dat[:,0]
                        sin_dat_imag = sin_dat[:,1]

                        sin_dat_comp = []
                        for i in range(len(sin_dat_real)):
                            sin_dat_comp.append(complex(sin_dat_real[i],sin_dat_imag[i]))

                        sin_dat_abs = []
                        for i in range(len(sin_dat_comp)):
                            sin_dat_abs.append(abs(sin_dat_comp[i]))


                        # Instrument State Dump in ASCII
                        sock.send("DSAS;" + NL)
                        time.sleep(15)
                        sock.send("++read eoi" + NL)
                        time.sleep(15)
                        sin_state = sock.recv(131072)
                        sin_state = [sin_state]
                        sin_state = sin_state[0].split('\r\n')
                        sin_state.pop()

                        # Retrieve useful information from the headers before generating plot
                        frq_strt = float(sin_state[94])
                        frq_end = float(sin_state[95])
                        frq_span = frq_end - frq_strt
                        frq_res = frq_span/801
                        x_cord = float(sin_header[11])
                        swp_rate = float(sin_state[68])
                        swp_time = float(sin_state[67])
                        swp_res = float(sin_state[69])
                        swp_intm = float(sin_state[70])

                        # Generate Frequency values based on start frequncy, frequncy span and frequency resolution
                        x_val_ss = np.arange(frq_strt,frq_strt+frq_span,frq_res)

                        session.run_type = request.vars.run_type
                        
                        if run_type == "PSD":
                            x_psd_pck = pickle.dumps(x_val)
                            x_psd_dat = pickle.dumps(dat)
                        else:
                            db.results.insert(time = strftime("%Y-%m-%d %H:%M:%S", localtime()), 
    										device = request.vars.device,
											run_type = request.vars.run_type,
                                            start_frequency = mSF_db_array[n],
                                            frequency_range = request.vars.frequency_range,
                                            source_level = request.vars.source_level,
                                            source_offset = request.vars.source_offset,
                                            sweep_rate = request.vars.sweep_rate)
                            x_pck = pickle.dumps(x_val_ss)
                            yr_pck = pickle.dumps(sin_dat_real)
                            ya_pck = pickle.dumps(sin_dat_abs)
                            db.data.insert(x = x_pck,
                                        yr = yr_pck,
                                        ya = ya_pck)


    elif form.errors:
        response.flash = 'Error, please try again'
    else:
        response.flash = 'Please fill in the data'

    rows = db(db.results.id==db.data.id).select(db.results.time,db.results.device,db.results.run_type,db.results.start_frequency,db.results.frequency_range,db.results.source_level,db.results.source_offset,db.results.sweep_rate,db.data.id)
    records=SQLTABLE(rows, headers='labels')
    return dict(form=form,vars=form.vars,records=records)

def plot_data():
    dataplot=db.data[request.args(0)] or redirect(URL(r=request,f='index'))
    x_val = dataplot.x
    y_abs = dataplot.ya
    y_real = dataplot.yr
    return image_mat(x_val,y_abs,y_real)

def image_mat(x,ya,yr):
    import sys
    import cStringIO
    from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
    from matplotlib.figure import Figure
    import matplotlib

    #run_type = request.vars.run_type

    #if (run_type == "PSD"):
    #    fig=Figure()
    #    ax=fig.add_subplot(111)
    #    x_val = session.x_val
    #    dat = session.dat
    #    ax.plot(x_val,dat)
    #    ax.set_xlabel("Hz")
    #    ax.set_ylabel("dB")
    #    ax.axis('tight')
    #    canvas=FigureCanvas(fig)
    #    response.headers['Content-Type']="image/png"
    #    stream=cStringIO.StringIO()
    #    canvas.print_png(stream)
    #elif (run_type == "Sweep"):
    x_val_ss = pickle.loads(x)
    sin_dat_abs = pickle.loads(ya)
    sin_dat_real = pickle.loads(yr)
    fig=Figure()
    ax=fig.add_subplot(211)
    ax.plot(x_val_ss,sin_dat_abs)
    ax.set_xlabel("Hz")
    ax.set_ylabel("Mag")
    ax.axis('tight')
    ay=fig.add_subplot(212)
    ay.plot(x_val_ss,sin_dat_real)
    ay.set_xlabel("Hz")
    ay.set_ylabel("Real")
    ay.axis('tight')
    sizer = request.function
    if sizer == 'index':
        rcParams['figure.figsize'] = 1.5,.75
    else:
        rcParams['figure.figsize'] = 10,8
    canvas=FigureCanvas(fig)
    response.headers['Content-Type']='image/png'
    stream=cStringIO.StringIO()
    canvas.print_png(stream)

    return stream.getvalue()

def plot_thumb():
    dataplot_th=db.data[request.args(0)] or redirect(URL(r=request,f='index'))
    x_val_th = dataplot_th.x
    y_abs_th = dataplot_th.ya
    y_real_th = dataplot_th.yr
    return image_mat_thumb(x_val_th,y_abs_th,y_real_th)

def image_mat_thumb(xth,yath,yrth):
    import sys
    import cStringIO
    from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
    from matplotlib.figure import Figure
    import matplotlib

    #run_type = request.vars.run_type

    #if (run_type == "PSD"):
    #    fig=Figure()
    #    ax=fig.add_subplot(111)
    #    x_val = session.x_val
    #    dat = session.dat
    #    ax.plot(x_val,dat)
    #    ax.set_xlabel("Hz")
    #    ax.set_ylabel("dB")
    #    ax.axis('tight')
    #    canvas=FigureCanvas(fig)
    #    response.headers['Content-Type']="image/png"
    #    stream=cStringIO.StringIO()
    #    canvas.print_png(stream)
    #elif (run_type == "Sweep"):
    x_val_ss = pickle.loads(xth)
    sin_dat_abs = pickle.loads(yath)
    sin_dat_real = pickle.loads(yrth)
    fig=Figure()
    ax=fig.add_subplot(211)
    ax.plot(x_val_ss,sin_dat_abs)
    ax.set_xlabel("Hz")
    ax.set_ylabel("Mag")
    ax.axis('tight')
    ay=fig.add_subplot(212)
    ay.plot(x_val_ss,sin_dat_real)
    ay.set_xlabel("Hz")
    ay.set_ylabel("Real")
    ay.axis('tight')
    sizer = request.function
    if sizer == 'index':
        rcParams['figure.figsize'] = 1.5,.75
    else:
        rcParams['figure.figsize'] = 10,8
    canvas=FigureCanvas(fig)
    response.headers['Content-Type']='image/png'
    stream=cStringIO.StringIO()
    canvas.print_png(stream)

    return stream.getvalue()

def calc_sweep_time():
    FRS_in = request.vars.frequency_range
    SW_in = request.vars.sweep_rate
    m_runs = request.vars.multiple_runs
    n_in = request.vars.number_runs

    n_runs = float(n_in)
    frq_range = float(FRS_in)
    swp_rate = float(SW_in)
    
    swp_t_single = round(((frq_range/swp_rate) + 90)/60, 1)
    
    if (m_runs == 'on'):
        swp_t_tot = swp_t_single * n_runs
    else:
        swp_t_tot = swp_t_single
        
    swp_t_s = str(swp_t_single)
    swp_t_t = str(swp_t_tot)

    return "jQuery('#run_time').html('Each run should take approximately %s minutes. The total sweep time is approximately %s minutes')" % (swp_t_s, swp_t_t)

def calc_sweep_func():
    SF_in = request.vars.frequency_start
    FRS_in = request.vars.frequency_range
    device = request.vars.device

    if device == "Medium":
        HOST = "172.21.7.68"
    else:
        HOST = "172.21.7.69"

    PORT = 1234
    NL = "\r\n"

    SF_pre = "SF "
    FRS_pre = "FRS "
    HZ_suf = " HZ;" + NL

    SF_var = binascii.a2b_qp(SF_pre+str(SF_in)+HZ_suf)
    FRS_var = binascii.a2b_qp(FRS_pre+str(FRS_in)+HZ_suf)

    sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM, socket.IPPROTO_TCP)
    sock.settimeout(0.5)
    sock.connect((HOST,PORT))

    sock.send("++mode 1" + NL)
    sock.send("++addr 20" + NL)
    sock.send("++auto 0" + NL)
    sock.send("++eoi 1" + NL)

    sock.send("AB;" + NL)
    time.sleep(0.5)
    sock.send("MSMD;" + NL)
    time.sleep(0.5)
    sock.send("SSIN;" + NL)
    time.sleep(0.5)
    sock.send("LNSW;" + NL)
    time.sleep(0.5)
    sock.send("MDSP;" + NL)
    time.sleep(0.5)
    sock.send("FRQR;" + NL)
    time.sleep(0.5)
    sock.send("SMES;" + NL)
    time.sleep(0.5)
    sock.send("FRSP;" + NL)
    time.sleep(0.5)
    sock.send("STTR;" + NL)
    time.sleep(0.5)
    sock.send("FREQ;" + NL)
    time.sleep(0.5)
    sock.send(SF_var)
    time.sleep(0.5)
    sock.send(FRS_var)
    time.sleep(0.5)
    sock.send("CORD;" + NL)
    time.sleep(0.5)
    sock.send("A;" + NL)
    time.sleep(0.5)
    sock.send("MAG;" + NL)
    time.sleep(0.5)
    sock.send("B;" + NL)
    time.sleep(0.5)
    sock.send("REAL;" + NL)
    time.sleep(0.5)
    sock.send("SCAL;" + NL)
    time.sleep(0.5)
    sock.send("AB;" + NL)
    time.sleep(0.5)
    sock.send("YASC;" + NL)
    time.sleep(0.5)
    sock.send("XASC;" + NL)
    time.sleep(0.5)
    sock.send("STTR;" + NL)
    time.sleep(0.5)
    sock.send("SWRT;" + NL)
    time.sleep(0.5)
    sock.send("SWRT?" + NL)
    time.sleep(6)
    sock.send("++read eoi" + NL)
    swp_el = sock.recv(64)
    swp_el = swp_el.split('\r\n')[0]

    return "jQuery('#gpib_sweep_rate').val(%s);\njQuery('#gpib_sweep_rate').blur();" % repr(swp_el)

def sweep_up_func():

    FRS_in = request.vars.frequency_range
    device = request.vars.device

    if device == "Medium":
        HOST = "172.21.7.68"
    else:
        HOST = "172.21.7.69"

    PORT = 1234
    NL = "\r\n"

    sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM, socket.IPPROTO_TCP)
    sock.settimeout(0.5)
    sock.connect((HOST,PORT))

    sock.send("SWRT;" + NL)
    time.sleep(0.5)
    sock.send("UP;" + NL)
    time.sleep(0.5)
    sock.send("SWRT?" + NL)
    time.sleep(6)
    sock.send("++read eoi" + NL)
    swp_el_up = sock.recv(64)
    swp_el_up = swp_el_up.split('\r\n')[0]

    return "jQuery('#gpib_sweep_rate').val(%s);\njQuery('#gpib_sweep_rate').blur();" % repr(swp_el_up)

def sweep_down_func():

    FRS_in = request.vars.frequency_range
    device = request.vars.device

    if device == "Medium":
        HOST = "172.21.7.68"
    else:
        HOST = "172.21.7.69"

    PORT = 1234
    NL = "\r\n"

    sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM, socket.IPPROTO_TCP)
    sock.settimeout(0.5)
    sock.connect((HOST,PORT))

    sock.send("SWRT;" + NL)
    time.sleep(0.5)
    sock.send("DOWN;" + NL)
    time.sleep(0.5)
    sock.send("SWRT?" + NL)
    time.sleep(6)
    sock.send("++read eoi" + NL)
    swp_el_down = sock.recv(64)
    swp_el_down = swp_el_down.split('\r\n')[0]

    return "jQuery('#gpib_sweep_rate').val(%s);\njQuery('#gpib_sweep_rate').blur();" % repr(swp_el_down)

def user():
    """
    exposes:
    http://..../[app]/default/user/login
    http://..../[app]/default/user/logout
    http://..../[app]/default/user/register
    http://..../[app]/default/user/profile
    http://..../[app]/default/user/retrieve_password
    http://..../[app]/default/user/change_password
    use @auth.requires_login()
        @auth.requires_membership('group name')
        @auth.requires_permission('read','table name',record_id)
    to decorate functions that need access control
    """
    return dict(form=auth())


def download():
    """
    allows downloading of uploaded files
    http://..../[app]/default/download/[filename]
    """
    return response.download(request,db)


def call():
    """
    exposes services. for example:
    http://..../[app]/default/call/jsonrpc
    decorate with @services.jsonrpc the functions to expose
    supports xml, json, xmlrpc, jsonrpc, amfrpc, rss, csv
    """
    session.forget()
    return service()



@auth.requires_signature()
def data():
    """
    http://..../[app]/default/data/tables
    http://..../[app]/default/data/create/[table]
    http://..../[app]/default/data/read/[table]/[id]
    http://..../[app]/default/data/update/[table]/[id]
    http://..../[app]/default/data/delete/[table]/[id]
    http://..../[app]/default/data/select/[table]
    http://..../[app]/default/data/search/[table]
    but URLs must be signed, i.e. linked with
      A('table',_href=URL('data/tables',user_signature=True))
    or with the signed load operator
      LOAD('default','data.load',args='tables',ajax=True,user_signature=True)
    """
    return dict(form=crud())
