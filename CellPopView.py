def CellPopView(adata, Cluster_Key, Cluster, Condition_Key, Conditions):
    #adata = annData object
    #Cluster_Key = key cluster labels are stored under in adata.obs
    #Cluster = specific cluster to be analyzed
    #Condition_Key = key the condition labels are stored under in adata.obs
    #Conditions = the two conditions to be compared, as a list
    
    #subset by cluster
    adata = adata[adata.obs[Cluster_Key].isin([Cluster])]
    
    #Split by Condition
    adata_1 = adata[adata.obs[Condition_Key].isin([Conditions[0]])]
    adata_2 = adata[adata.obs[Condition_Key].isin([Conditions[1]])]
    
    #Extract the data
    Table_1 = adata_1.to_df()
    Table_2 = adata_2.to_df()
    
    Table_1 = Table_1.transpose() #now genes = rows, cells = columns
    Expression_1 = Table_1.mean(axis=1) #extracts the average of every row, e.g. average expression of every gene

    Table_2 = Table_2.transpose()
    Expression_2 = Table_2.mean(axis=1) 
    
    #plot and return graph
    plt.scatter(Expression_1,Expression_2, label = "stars", color = "green", 
                marker = "*",  s =30) 
    
    plt.title('Cluster ' + str(Cluster)) 
    
    plt.xlabel(Conditions[0])
    plt.ylabel(Conditions[1], rotation = 0)
    
    CellPopPlot = plt.show()
    
    return(CellPopPlot)
