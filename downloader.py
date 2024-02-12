"""
This is the downloader for the repository. There are only the raw data scripts in here now but in the future I'll add other scripts
If it is necessary to get other scripts, just add the raw path to the list of scripts. 
"""
import requests

if __name__ == '__main__':
    base_path = 'https://raw.githubusercontent.com/petegombert/cloud_counter/main/'
    script_list = ['oned_corr','get_data','mast_plot','go_stats','build_ann_anom','coarsify','diveinnewdata','jja_global_avg','overview','indexes','bam_analysis']

    for script in script_list:
        URL = base_path+script+'.py'
        with requests.Session() as session:
            FILENAME = script+'.py'
            r1 = session.request('get',URL)
            r = session.get(r1.url)
            r.raise_for_status()
          
            f = open(FILENAME,'wb')
            f.write(r.content)
            f.close()
            print('Downloaded: '+script+'.py')
