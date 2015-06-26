from BaseSpacePy.api.BaseSpaceAPI import BaseSpaceAPI
import argparse,sys,os

def ParseArg():
    p=argparse.ArgumentParser( description="")
    p.add_argument('-k','--key',type=str,default="4f7366779990451799fc491d4f6f51b5",help="the client_key, default: 4f7366779990451799fc491d4f6f51b5")
    p.add_argument('-s','--secret',type=str,default="c88e8c58e5814d04b3c662dc199693ab",help="the client_serect, default: c88e8c58e5814d04b3c662dc199693ab")
    p.add_argument('-t','--token',type=str,default="639f9af2f031415cb89bfeef18716b72",help="the accessToken for the app [NECESSARY], default: 639f9af2f031415cb89bfeef18716b72")
    p.add_argument('-p','--project',type=str,required=True,help="query project full name.")
    p.add_argument('-d','--directory',type=str,default='./',help="download file folder, default: current folder")
    if len(sys.argv)==1:
        print >>sys.stderr,p.print_help()
        exit(0)
    return p.parse_args()

def Main():
    args = ParseArg()
    client_key = args.key
    client_secret = args.secret
    token = args.token
    folder = args.directory
    BaseSpaceUrl = 'https://api.basespace.illumina.com/'
    myAPI = BaseSpaceAPI(client_key,client_secret,BaseSpaceUrl,"v1pre3","AA",token)
    user = myAPI.getUserById('current')
    print >> sys.stderr, "\nUser name: %s\n"%(str(user))
    Projects = myAPI.getProjectByUser()
    Found = False
    for p in Projects:
        if p.Name == args.project:
            print >>sys.stderr, "  Find project %s with ID: %s. "%(p.Name, p.Id)
            Project = p
            Found = True
            break
    if not Found: 
        print >>sys.stderr, "  Could not find project %s, from user %s, please check your token." %(args.project,str(user))
        sys.exit(0)

    Samples=Project.getSamples(myAPI)
    print >> sys.stderr, "Samples for this project: " + str(Samples)
    for s in Samples:
        print >>sys.stderr,"  Downloading files in sample " + str(s)
        subfolder=folder+"/"+str(s)
        if not os.path.exists(subfolder):
            os.makedirs(subfolder)
        for f in s.getFiles(myAPI):
            print >> sys.stderr,"    "+str(f)
            f.downloadFile(myAPI,subfolder)
    
if __name__=="__main__":
    Main()
