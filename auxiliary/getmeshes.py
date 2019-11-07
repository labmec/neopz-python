import requests

def Simple2DMesh():
    url = 'https://raw.githubusercontent.com/labmec/neopz-python/master/tests/geometric-mesh/simple_2D_coarse.msh'
    r = requests.get(url, allow_redirects=True)
    open('simple_2D_coarse.msh', 'wb').write(r.content)
    return 'simple_2D_coarse.msh'

def DownloadFromURL(url):
    spliturl = url.split('/')
    name = spliturl[len(spliturl)-1]
    splitname = name.split('.')
    extension = splitname[len(splitname)-1]
    if extension != 'msh':
        print("The file must be .msh extension!")
    else:
        r = requests.get(url, allow_redirects=True)
        open(name, 'wb').write(r.content)
        return name


