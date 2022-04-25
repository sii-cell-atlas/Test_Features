def CellPopView(data):
    
    adata = createData(data)

    #log_normalize
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    #subset by cluster
    Cluster_Key = data['ClusterKey']
    
    Cluster = data['Cluster']
    
    adata = adata[adata.obs[Cluster_Key].isin([Cluster])]

    #Split by Condition
    Condition_Key = data['ConditionKey']
    
    Condition1 = data['Cond1']
    
    Condition2 = data['Cond2']
    
    adata_1 = adata[adata.obs[Condition_Key].isin([Condition1])]
    adata_2 = adata[adata.obs[Condition_Key].isin([Condition2])]
    
    
    #Extract the data
    Table_1 = adata_1.to_df()
  
    Table_2 = adata_2.to_df()
    
    Table_1 = Table_1.transpose() #now genes = rows, cells = columns
    Expression_1 = Table_1.mean(axis=1) #extracts the average of every row, e.g. average expression of every gene
    
    Table_2 = Table_2.transpose()
    Expression_2 = Table_2.mean(axis=1) #average expression of every gene in Cluster 3, Condition S
    
    plt.scatter(Expression_1,Expression_2, label = "stars", color = "black", 
            marker = ".",  s =30) 
    
    plt.title(Cluster_Key + ": " + Cluster) 
    
    plt.grid()

    plt.xlabel(Condition1)
    plt.ylabel(Condition2)
    
    CellPopPlot = plt.gcf()
    
    return iostreamFig(CellPopPlot)
