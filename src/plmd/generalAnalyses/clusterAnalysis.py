#!/usr/bin/python
import plmd.plotData as myPlot

class clusterBase():
    
    # Constructor
    def __init__( self , dataFile ):
        
        # Get the cluster vs. frame information
        self.clusterInfo = []
        with open( dataFile , "r") as fi:
            next(fi)
            for line in fi:
                self.clusterInfo.append( line.split()[1] )
        
    # A function for separating a dataset based on cluster file
    def separateDataSet( self, dataFile , outFiles , xColumn = 0 ):
        data = {}
        clusterLen = len(self.clusterInfo)
        clusters = []
        with open(dataFile, "r") as fi:
            next(fi)
            for line in fi:
                temp = line.split()
                
                # Get the cluster
                if clusterLen >= int(temp[0]):
                    
                    # Clusters
                    clusterNr = self.clusterInfo[ int(temp[0])-1 ]
                    if clusterNr not in clusters:
                        clusters.append(clusterNr)           
                    
                    # Add this cluster for each component in data file
                    for i in range(1, len(temp)):
                        
                        # The identifier for this dataset
                        cluster = "d"+str(i)+"_c"+clusterNr
                        
                        # Check if it already exists
                        if cluster not in data:
                            data[ cluster ] = { "x":[], "y":[] }
                        
                        # Add entry
                        data[ cluster ]["x"].append( temp[ xColumn ] )
                        data[ cluster ]["y"].append( temp[ i ] )
    
        print "==================="
        print clusters
        # Save it in appropriate files
        for cluster, dataset in data.iteritems():
            with open( outFiles+str(cluster) , "w") as fo:
                for i in range( 0, len(dataset["x"]) ):
                    fo.write( dataset["x"][i]+"\t"+dataset["y"][i]+"\n" )
          
        return len(clusters)

