#Similarity between two point is measured as euclidean distance(L2 Norm)
get_euclidean_distance = function(data_point1, data_point_2) {
  #print(dim(data_point1))
  total_dist = 0
  no_of_features = ncol(data_point1)
  if(no_of_features != ncol(data_point_2)) {
    return (total_dist)
  }
  else {
    for(i in 1:no_of_features) {
      total_dist = total_dist + (data_point1[i] - data_point_2[i])^2
    }
    total_dist = total_dist^(0.5)
  }
  #print(total_dist)
  return (total_dist)
}


#original Proximity(Similarity Matrix of the data points)
calculate_distance = function(datapoints) {
  
  similarity_matrix = matrix(nrow=nrow(datapoints), ncol=nrow(datapoints), byrow = TRUE)
  #print(dim(similarity_matrix))

  
  for(i in 1:nrow(datapoints)) {
    for(j in 1:i) {
      if (i == j){
        similarity_matrix[i,j] = NA
      } else {
        x = as.matrix(datapoints[i,])
        y = as.matrix(datapoints[j,])
        #print(y)
        euc_distance = as.numeric(get_euclidean_distance(x, y))
        similarity_matrix[i,j] = euc_distance
        similarity_matrix[j,i] = euc_distance
      }
    }  
  }
  #print(similarity_matrix)
  return(similarity_matrix)
}

# Single(min) linkage method
single_method = function(distance_matrix, K, cluster) {
  
  no_row = nrow(distance_matrix)
  valid_data_points = rep(1,no_row)  #intially all points to be considered and thus valid
  
  #Max Rows of merge matrix
  temp_cluster = matrix(NA, no_row*2, no_row*2)
  temp_cluster[1:nrow(distance_matrix),1:ncol(distance_matrix)] = distance_matrix
  
  distance_matrix = temp_cluster
  
  for (idx in 1:(no_row-K)) {
    new_row = no_row + idx

    join_cluster = which(distance_matrix==min(distance_matrix, na.rm=TRUE), arr.ind=TRUE)
    merge_x = join_cluster[nrow(join_cluster),1]
    merge_y = join_cluster[nrow(join_cluster),2]
    cluster[[new_row]] = c(cluster[[merge_x]],cluster[[merge_y]])  
    
    new_cluster_pt = cluster[[length(cluster)]]
    for (pt in c(1:(length(cluster)-1))[-new_cluster_pt]) {
      merge_cluster = c()
      for (j in cluster[[pt]]) {
        for (k in new_cluster_pt) {
          merge_cluster = append(merge_cluster, temp_cluster[j,k])
        }
      }
      temp_cluster[length(cluster),pt] =  temp_cluster[pt,length(cluster)] = min(merge_cluster)
    }
    
    valid_data_points[join_cluster[nrow(join_cluster),]] = 0
    valid_data_points[new_row] = 1
    
    distance_matrix[,join_cluster[nrow(join_cluster),]] = distance_matrix[join_cluster[nrow(join_cluster),],] = NA
    distance_matrix[(new_row),1:(new_row)] = temp_cluster[(new_row),1:(new_row)]
    distance_matrix[1:(new_row),(new_row)] = temp_cluster[1:(new_row),(new_row)]
  }
  return (cluster[valid_data_points==1])
}


# Complete(max) linkage method
complete_method = function(distance_matrix, K, cluster) {
  no_row = nrow(distance_matrix)
  valid_data_points = c(rep(1,no_row))  #intially all points to be considered and thus valid
  
  #Max Rows of merge matrix
  temp_cluster = matrix(NA, no_row*2, no_row*2)
  temp_cluster[1:nrow(distance_matrix),1:ncol(distance_matrix)] = distance_matrix
  
  distance_matrix = temp_cluster
  
  for (idx in 1:(no_row-K)) {
    new_row = no_row + idx
    
    join_cluster = which(distance_matrix==min(distance_matrix, na.rm=TRUE), arr.ind=TRUE)
    merge_x = join_cluster[nrow(join_cluster),1]
    merge_y = join_cluster[nrow(join_cluster),2]
    cluster[[new_row]] = c(cluster[[merge_x]],cluster[[merge_y]])  
    
    new_cluster_pt = cluster[[length(cluster)]]
    for (pt in c(1:(length(cluster)-1))[-new_cluster_pt]) {
      merge_cluster = c()
      for (j in cluster[[pt]]) {
        for (k in new_cluster_pt) {
          merge_cluster = append(merge_cluster, temp_cluster[j,k])
        }
      }
      temp_cluster[length(cluster),pt] =  temp_cluster[pt,length(cluster)] = max(merge_cluster)
    }
    
    valid_data_points[join_cluster[nrow(join_cluster),]] = 0
    valid_data_points[new_row] = 1
    
    distance_matrix[,join_cluster[nrow(join_cluster),]] = distance_matrix[join_cluster[nrow(join_cluster),],] = NA
    distance_matrix[(new_row),1:(new_row)] = temp_cluster[(new_row),1:(new_row)]
    distance_matrix[1:(new_row),(new_row)] = temp_cluster[1:(new_row),(new_row)]
  }
  return (cluster[valid_data_points==1])
} 


# Group Average linkage method
avg_method = function(distance_matrix, K, cluster) {
  no_row = nrow(distance_matrix)
  valid_data_points = c(rep(1,no_row))  #intially all points to be considered and thus valid
  
  #Max Rows of merge matrix
  temp_cluster = matrix(NA, no_row*2, no_row*2)
  temp_cluster[1:nrow(distance_matrix),1:ncol(distance_matrix)] = distance_matrix
  
  distance_matrix = temp_cluster
  
  for (idx in 1:(no_row-K)) {
    new_row = no_row + idx
    
    join_cluster = which(distance_matrix==min(distance_matrix, na.rm=TRUE), arr.ind=TRUE)
    merge_x = join_cluster[nrow(join_cluster),1]
    merge_y = join_cluster[nrow(join_cluster),2]
    cluster[[new_row]] = c(cluster[[merge_x]],cluster[[merge_y]])  
    
    new_cluster_pt = cluster[[length(cluster)]]
    for (pt in c(1:(length(cluster)-1))[-new_cluster_pt]) {
      merge_cluster = c()
      for (j in cluster[[pt]]) {
        for (k in new_cluster_pt) {
          merge_cluster = append(merge_cluster, temp_cluster[j,k])
        }
      }
      temp_cluster[length(cluster),pt] =  temp_cluster[pt,length(cluster)] = sum(merge_cluster)/(length(cluster[[pt]])*length(new_cluster_pt))
    }
    
    valid_data_points[join_cluster[nrow(join_cluster),]] = 0
    valid_data_points[new_row] = 1
    
    distance_matrix[,join_cluster[nrow(join_cluster),]] = distance_matrix[join_cluster[nrow(join_cluster),],] = NA
    distance_matrix[(new_row),1:(new_row)] = temp_cluster[(new_row),1:(new_row)]
    distance_matrix[1:(new_row),(new_row)] = temp_cluster[1:(new_row),(new_row)]
  }
  return (cluster[valid_data_points==1])
} 

# Centroid Linkage method
centroid_method = function(data_set, distance_matrix, K, cluster) {
  no_row = nrow(distance_matrix)
  valid_data_points = c(rep(1,no_row))  #intially all points to be considered and thus valid
  
  #Max Rows of merge matrix
  temp_cluster = matrix(NA, no_row*2, no_row*2)
  temp_cluster[1:nrow(distance_matrix),1:ncol(distance_matrix)] = distance_matrix
  
  distance_matrix = temp_cluster
  
  for (idx in 1:(no_row-K)) {
    new_row = no_row + idx
    
    join_cluster = which(distance_matrix==min(distance_matrix, na.rm=TRUE), arr.ind=TRUE)
    merge_x = join_cluster[nrow(join_cluster),1]
    merge_y = join_cluster[nrow(join_cluster),2]
    cluster[[new_row]] = c(cluster[[merge_x]],cluster[[merge_y]])  
    
    new_cluster_pt = cluster[[length(cluster)]]
    for (pt in c(1:(length(cluster)-1))[-new_cluster_pt]) {
      merge_cluster = c()
      for (j in cluster[[pt]]) {
        for (k in new_cluster_pt) {
          merge_cluster = append(merge_cluster, temp_cluster[j,k])
        }
      }
      centriod_dist = abs(mean(merge_cluster))
      temp_cluster[length(cluster),pt] =  temp_cluster[pt,length(cluster)] = centriod_dist
    }
    
    valid_data_points[join_cluster[nrow(join_cluster),]] = 0
    valid_data_points[new_row] = 1
    
    distance_matrix[,join_cluster[nrow(join_cluster),]] = distance_matrix[join_cluster[nrow(join_cluster),],] = NA
    distance_matrix[(new_row),1:(new_row)] = temp_cluster[(new_row),1:(new_row)]
    distance_matrix[1:(new_row),(new_row)] = temp_cluster[1:(new_row),(new_row)]
  }
  return (cluster[valid_data_points==1])
}
