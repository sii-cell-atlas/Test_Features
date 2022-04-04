def CellPopTable(data):
  adata = createData(data)

  #subset by cluster
  Cluster_Key = data['ClusterKey']
    
  Cluster = data['Cluster']
    
  adata = adata[adata.obs[Cluster_Key].isin([Cluster])]

  #remove extraneous conditions
  
  Condition_Key = data['ConditionKey']
    
  Condition1 = data['Cond1']
    
  Condition2 = data['Cond2']
    
  adata = adata[adata.obs[Condition_Key].isin([Condition1,Condition2])]

  #Differential Expression Analysis
  res = de.test.t_test(adata,grouping=Condition_Key)

  deg = res.summary()
  deg = deg.sort_values(by=['qval']).loc[:,['gene','log2fc','pval','qval']]
  deg['log2fc'] = -1 * deg['log2fc']

  gInfo = getVar(data)
  deg.index = deg['gene']
  deg = pd.concat([deg,gInfo],axis=1,sort=False)

  restable = deg.to_html()
  
  ppr.pprint("Basic HTML Table created")

  return json.dumps(restable)
