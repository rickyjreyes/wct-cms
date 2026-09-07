"""Retrieve all four public ROOT files and verify published Adler-32 values."""
import pathlib,json,subprocess,concurrent.futures,zlib

def main():
    root=pathlib.Path('data/audit_inputs');root.mkdir(parents=True,exist_ok=True);jobs=[]
    for record,prefix in [('30555',('127C','183B')),('30522',('05DD','209D'))]:
        for f in json.load(open(f'audit/2026-09-06/baseline/cms{record}.json'))['metadata']['_file_indices'][0]['files']:
            if f['filename'].startswith(prefix):jobs.append(f)
    def fetch(f):
        path=root/f['filename'];url='https://opendata.cern.ch'+f['uri'].split('eospublic.cern.ch/')[1]
        if not path.exists() or path.stat().st_size!=f['size']:
            subprocess.run(['curl','-L','--fail','--retry','3','-C','-',url,'-o',str(path)],check=True)
        checksum=1
        with path.open('rb') as stream:
            while chunk:=stream.read(8*1024**2):checksum=zlib.adler32(chunk,checksum)
        actual=f'adler32:{checksum:08x}'
        assert actual==f['checksum'],(path,actual,f['checksum'])
        return dict(path=str(path),url=url,bytes=path.stat().st_size,checksum=actual)
    with concurrent.futures.ThreadPoolExecutor(4) as pool:results=list(pool.map(fetch,jobs))
    pathlib.Path('audit/2026-09-06/download_checksums.json').write_text(json.dumps(results,indent=2))
if __name__=='__main__':main()
