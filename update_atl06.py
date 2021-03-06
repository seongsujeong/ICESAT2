#!/usr/bin/env python
#
#  Functionality:
#   - Searches ATL06 data from 
# Run this script without input argument to see the usage instructions and examples.
#
# Code by Seongsu Jeong at UC Irvine
#
# This code is under MIT license, same as the original code by NSIDC. More information of the original code is as below:
# ------------------------------------------------------------------------------
# NSIDC Data Download Script
# Tested in Python 2.7 and Python 3.4, 3.6, 3.7
#
# To run the script at a Linux, macOS, or Cygwin command-line terminal:
#   $ python nsidc-data-download.py
#
# On Windows, open Start menu -> Run and type cmd. Then type:
#     python nsidc-data-download.py
#
# The script will first search Earthdata for all matching files.
# You will then be prompted for your Earthdata username/password
# and the script will download the matching files.
#
# If you wish, you may store your Earthdata username/password in a .netrc
# file in your $HOME directory and the script will automatically attempt to
# read this file. The .netrc file should have the following format:
#    machine urs.earthdata.nasa.gov login myusername password mypassword
# where 'myusername' and 'mypassword' are your Earthdata credentials.

from __future__ import print_function

import base64
import itertools
import json
import netrc
import ssl
import sys
import os
from getpass import getpass
import datetime
import multiprocessing

try:
    from urllib.parse import urlparse
    from urllib.request import urlopen, Request, build_opener, HTTPCookieProcessor
    from urllib.error import HTTPError, URLError
except ImportError:
    from urlparse import urlparse
    from urllib2 import urlopen, Request, HTTPError, URLError, build_opener, HTTPCookieProcessor

#short_name = 'ATL06'
#version = '002'
#time_start = '2018-10-14T00:00:00Z'
#time_end = '2020-05-05T12:45:02Z'
#bounding_box = ''
#polygon = '-66.15691061690228,-62.59317615574119,-78.35830038826326,-69.16349415975444,-98.11823949593786,-69.24053171586348,-114.38250197776,-70.93841790150154,-141.3003778231269,-72.04893903612022,-165.49384821924374,-75.33106539393404,177.5504098292959,-75.38247296471017,173.3497037884124,-69.10725156162574,156.76631499780524,-66.54787716726221,141.41240398973923,-63.53516625896447,130.42006340408514,-62.548492593715814,115.00200733067035,-62.88621697482232,97.46917734176803,-61.72560776837533,80.82440771666971,-63.97358742480969,75.27827929633031,-66.67270751538402,56.808013940748275,-62.05277315524081,10.915370276577704,-67.35123764014733,-8.454068426794686,-68.0764961239038,-28.382829506614637,-71.54876687779381,-40.88764485852986,-74.8115704633071,-53.9855375373903,-70.0638400601326,-49.967789567167856,-62.05158641781778,-55.06487847988186,-59.58983556733419,-61.055800565890216,-60.946797287081665,-66.15691061690228,-62.59317615574119'
#filename_filter = ''
#url_list = []


#pre-defined polygon
dict_polygon={
    'ant':'-66.15691061690228,-62.59317615574119,-78.35830038826326,-69.16349415975444,-98.11823949593786,-69.24053171586348,-114.38250197776,-70.93841790150154,-141.3003778231269,-72.04893903612022,-165.49384821924374,-75.33106539393404,177.5504098292959,-75.38247296471017,173.3497037884124,-69.10725156162574,156.76631499780524,-66.54787716726221,141.41240398973923,-63.53516625896447,130.42006340408514,-62.548492593715814,115.00200733067035,-62.88621697482232,97.46917734176803,-61.72560776837533,80.82440771666971,-63.97358742480969,75.27827929633031,-66.67270751538402,56.808013940748275,-62.05277315524081,10.915370276577704,-67.35123764014733,-8.454068426794686,-68.0764961239038,-28.382829506614637,-71.54876687779381,-40.88764485852986,-74.8115704633071,-53.9855375373903,-70.0638400601326,-49.967789567167856,-62.05158641781778,-55.06487847988186,-59.58983556733419,-61.055800565890216,-60.946797287081665,-66.15691061690228,-62.59317615574119',
    'gre':'-63.8175205605513,74.85800169902663,-58.10401026258381,74.11346163350383,-55.58451598473744,68.70641419988947,-55.30914889055835,66.07334639497095,-49.143717958285436,60.29373468235471,-42.44163050594869,58.49876285849441,-37.90451433875775,64.25833842926329,-32.70655689680123,65.56887711515871,-30.40240201554316,67.25749853151709,-19.901153018000837,69.32341602889063,-19.438816813828012,72.88795386452698,-13.756285225494835,76.51097776550014,-14.304251539831792,79.10123203372567,-6.835693048982446,82.62890177411286,-36.78298430582492,84.49174062865558,-72.44103785337352,81.15930332474665,-76.4033224631926,77.65054330771864,-72.92635571028501,75.53960803389364,-67.57821420721443,75.11326309635521,-63.8175205605513,74.85800169902663'
}

CMR_URL = 'https://cmr.earthdata.nasa.gov'
URS_URL = 'https://urs.earthdata.nasa.gov'
CMR_PAGE_SIZE = 2000
CMR_FILE_URL = ('{0}/search/granules.json?provider=NSIDC_ECS'
                '&sort_key[]=start_date&sort_key[]=producer_granule_id'
                '&scroll=true&page_size={1}'.format(CMR_URL, CMR_PAGE_SIZE))


def get_username():
    username = ''

    # For Python 2/3 compatibility:
    try:
        do_input = raw_input  # noqa
    except NameError:
        do_input = input

    while not username:
        try:
            username = do_input('Earthdata username: ')
        except KeyboardInterrupt:
            quit()
    return username


def get_password():
    password = ''
    while not password:
        try:
            password = getpass('password: ')
        except KeyboardInterrupt:
            quit()
    return password


def get_credentials(url):
    """Get user credentials from .netrc or prompt for input."""
    credentials = None
    errprefix = ''
    try:
        info = netrc.netrc()
        username, account, password = info.authenticators(urlparse(URS_URL).hostname)
        errprefix = 'netrc error: '
    except Exception as e:
        if (not ('No such file' in str(e))):
            print('netrc error: {0}'.format(str(e)))
        username = None
        password = None

    while not credentials:
        if not username:
            username = get_username()
            password = get_password()
        credentials = '{0}:{1}'.format(username, password)
        credentials = base64.b64encode(credentials.encode('ascii')).decode('ascii')

        if url:
            try:
                req = Request(url)
                req.add_header('Authorization', 'Basic {0}'.format(credentials))
                opener = build_opener(HTTPCookieProcessor())
                opener.open(req)
            except HTTPError:
                print(errprefix + 'Incorrect username or password')
                errprefix = ''
                credentials = None
                username = None
                password = None

    return credentials


def build_version_query_params(version):
    desired_pad_length = 3
    if len(version) > desired_pad_length:
        print('Version string too long: "{0}"'.format(version))
        quit()

    version = str(int(version))  # Strip off any leading zeros
    query_params = ''

    while len(version) <= desired_pad_length:
        padded_version = version.zfill(desired_pad_length)
        query_params += '&version={0}'.format(padded_version)
        desired_pad_length -= 1
    return query_params


def build_cmr_query_url(short_name, version, time_start, time_end,
                        bounding_box=None, polygon=None,
                        filename_filter=None):
    params = '&short_name={0}'.format(short_name)
    params += build_version_query_params(version)
    params += '&temporal[]={0},{1}'.format(time_start, time_end)
    if polygon:
        params += '&polygon={0}'.format(polygon)
    elif bounding_box:
        params += '&bounding_box={0}'.format(bounding_box)
    if filename_filter:
        option = '&options[producer_granule_id][pattern]=true'
        params += '&producer_granule_id[]={0}{1}'.format(filename_filter, option)
    return CMR_FILE_URL + params


def cmr_download(urls,path=None,overwrite=False):
    """Download files from list of urls."""

    if path==None:
        path_out=os.getcwd()
    else:
        path_out=path
    if not urls:
        return

    url_count = len(urls)
    print('Downloading {0} files...'.format(url_count))
    credentials = None

    for index, url in enumerate(urls, start=1):
        if not credentials and urlparse(url).scheme == 'https':
            credentials = get_credentials(url)

        filename = '{}/{}'.format(path_out,url.split('/')[-1])
        if os.path.exists(filename) and (not overwrite):
            print('File exists:',filename)
            continue
        else:
            print('{0}/{1}: {2}'.format(str(index).zfill(len(str(url_count))),
                                        url_count,
                                        filename))

            try:
                # In Python 3 we could eliminate the opener and just do 2 lines:
                # resp = requests.get(url, auth=(username, password))
                # open(filename, 'wb').write(resp.content)
                req = Request(url)
                if credentials:
                    req.add_header('Authorization', 'Basic {0}'.format(credentials))
                opener = build_opener(HTTPCookieProcessor())
                data = opener.open(req).read()
                open(filename, 'wb').write(data)
            except HTTPError as e:
                print('HTTP error {0}, {1}'.format(e.code, e.reason))
            except URLError as e:
                print('URL error: {0}'.format(e.reason))
            except IOError:
                raise
            except KeyboardInterrupt:
                quit()

def cmr_download_parallel(urls,numworker,path=None,overwrite=False):
    #TODO: TEST THIS MORE. Make it work.
    #slice urls
    arg_seg=[None]*numworker
    for i in range(numworker):
        arg_seg[i]=(urls[i::3],path,overwrite)
    
    print(arg_seg[0])
    print(arg_seg[1])
    pin=multiprocessing.Pool(numworker)
    #pin.map(self.Gaussian_interp_single_rev2,list_cubeid)
    pin.starmap(cmr_download,arg_seg)
    


def cmr_filter_urls(search_results):
    """Select only the desired data files from CMR response."""
    if 'feed' not in search_results or 'entry' not in search_results['feed']:
        return []

    entries = [e['links']
               for e in search_results['feed']['entry']
               if 'links' in e]
    # Flatten "entries" to a simple list of links
    links = list(itertools.chain(*entries))

    urls = []
    unique_filenames = set()
    for link in links:
        if 'href' not in link:
            # Exclude links with nothing to download
            continue
        if 'inherited' in link and link['inherited'] is True:
            # Why are we excluding these links?
            continue
        if 'rel' in link and 'data#' not in link['rel']:
            # Exclude links which are not classified by CMR as "data" or "metadata"
            continue

        if 'title' in link and 'opendap' in link['title'].lower():
            # Exclude OPeNDAP links--they are responsible for many duplicates
            # This is a hack; when the metadata is updated to properly identify
            # non-datapool links, we should be able to do this in a non-hack way
            continue

        filename = link['href'].split('/')[-1]
        if filename in unique_filenames:
            # Exclude links with duplicate filenames (they would overwrite)
            continue
        unique_filenames.add(filename)

        urls.append(link['href'])

    return urls


def cmr_search(short_name, version, time_start, time_end,
               bounding_box='', polygon='', filename_filter=''):
    """Perform a scrolling CMR query for files matching input criteria."""
    cmr_query_url = build_cmr_query_url(short_name=short_name, version=version,
                                        time_start=time_start, time_end=time_end,
                                        bounding_box=bounding_box,
                                        polygon=polygon, filename_filter=filename_filter)
    print('Querying for data:\n\t{0}\n'.format(cmr_query_url))

    cmr_scroll_id = None
    ctx = ssl.create_default_context()
    ctx.check_hostname = False
    ctx.verify_mode = ssl.CERT_NONE

    try:
        urls = []
        while True:
            req = Request(cmr_query_url)
            if cmr_scroll_id:
                req.add_header('cmr-scroll-id', cmr_scroll_id)
            response = urlopen(req, context=ctx)
            if not cmr_scroll_id:
                # Python 2 and 3 have different case for the http headers
                headers = {k.lower(): v for k, v in dict(response.info()).items()}
                cmr_scroll_id = headers['cmr-scroll-id']
                hits = int(headers['cmr-hits'])
                if hits > 0:
                    print('Found {0} matches.'.format(hits))
                else:
                    print('Found no matches.')
            search_page = response.read()
            search_page = json.loads(search_page.decode('utf-8'))
            url_scroll_results = cmr_filter_urls(search_page)
            if not url_scroll_results:
                break
            if hits > CMR_PAGE_SIZE:
                print('.', end='')
                sys.stdout.flush()
            urls += url_scroll_results

        if hits > CMR_PAGE_SIZE:
            print()
        return urls
    except KeyboardInterrupt:
        quit()


def main(short_name, version, time_start, time_end, bounding_box='',polygon=None, filename_filter=None, url_list=None):
    #global short_name, version, time_start, time_end, bounding_box, \
    #    polygon, filename_filter, url_list

    # Supply some default search parameters, just for testing purposes.
    # These are only used if the parameters aren't filled in up above.
    
    url_list = cmr_search(short_name, version,
                              time_start, time_end,
                              bounding_box=bounding_box,
                              polygon=polygon, filename_filter=filename_filter)

    cmr_download(url_list)


if __name__ == '__main__':
    str_usage='''
    Usage example:
    update_atl06.py -r ant -t0 2001-01-01 -t1 2020-05-05 -o /downloads'
    update_atl06.py -r gre -t0 2001-01-01 -t1 2020-05-05 -o /downloads'
    update_atl06.py -r gre -w 2 -t0 2001-01-01 -t1 2020-05-05 -o /downloads'
    update_atl06.py -r ant -t0 2001-01-01 -t1 2020-05-05 -o /downloads --overwrite'
    update_atl06.py -r ant -v 003 -t0 2001-01-01 -t1 2020-05-05 -o /downloads'

    Options:
    -r: Region [ant|gre] default: ant
    -v: Data version, default: 003
    -w: Number of workers i.e. concurrent download streams. Default:1
    -t0: Starting date, default: 2018-01-01
    -t1: Ending time, default: TODAY
    -o: path to save the .h5 and xml files. Default: current working directoy
    
    --overwrite: Use this option to overwrite the files

    Based on NSIDC ICESAT-2 data retrieval script
    Modified by Seongsu Jeong at University of California, Irvine

    NOTE about parallel mode
    - Parallel mode is currently under experiment.
    - Please be a good user to the NASA DAAC. DO NOT spawn too many download worker as it might cause the server to overload.

    '''
    #default paramaters
    str_shortname='ATL06'
    str_version='003'
    str_boundingbox=''
    str_filename_filter=''
    numworker=1

    str_region='ant'
    time_start='2018-10-14T00:00:00Z'
    time_end=datetime.datetime.now().strftime('%Y-%m-%dT%H:%M:%SZ')
    path_out=os.getcwd()
    flag_overwrite=False
    
    #parse the input arguments
    argid=1
    
    if len(sys.argv)==1:
        print(str_usage)
    else:
        while argid<len(sys.argv):
            if sys.argv[argid]=='-r':
                str_region=sys.argv[argid+1]
                argid+=2
            elif sys.argv[argid]=='-v':
                str_version=sys.argv[argid+1]
                argid+=2
            elif sys.argv[argid]=='-w':
                numworker=int(sys.argv[argid+1])
                argid+=2
            elif sys.argv[argid]=='-t0':
                time_start='{}T00:00:00Z'.format(sys.argv[argid+1])
                argid+=2
            elif sys.argv[argid]=='-t1':
                time_end='{}T23:59:59Z'.format(sys.argv[argid+1])
                argid+=2
            elif sys.argv[argid]=='-o':
                path_out=sys.argv[argid+1]
                argid+=2
            elif sys.argv[argid]=='--overwrite':
                flag_overwrite=True
                argid+=2
            else:
                print('Cannot understand the input argument:',sys.argv[argid])
                exit(1)

        #main()

            #determine the polygon
        str_polygon=dict_polygon[str_region]


        url_list = cmr_search(str_shortname, str_version,
                            time_start, time_end,
                            bounding_box=str_boundingbox,
                            polygon=str_polygon, filename_filter=str_filename_filter)

        if numworker==1:
            cmr_download(url_list,path=path_out,overwrite=flag_overwrite)
        else:
            print('numworker:',numworker)
            cmr_download_parallel(url_list,numworker,path=path_out,overwrite=flag_overwrite)
