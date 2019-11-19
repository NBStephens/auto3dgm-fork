write_aligned_files <-
function(ds, ga, input_dir, output_dir, scale){
  #rotate the individual raw files.  Use the rotation matrices from global alignment.  
  meshes.flipped<-c()
  for (kk in 1: ds$n){
    mesh_fn = file.path(input_dir, ds$ids[kk])
    mesh_fn_type = substring(tolower(mesh_fn), nchar(mesh_fn)-3, nchar(mesh_fn))
    
    if(mesh_fn_type == ".off"){
      VF = read_off(mesh_fn)
    }else if(mesh_fn_type == ".ply"){
      VF = read_ply(mesh_fn)
    }else{
      stop("Error. Unrecognized mesh file type.")
    }
    
    #VF = read_off(file.path(input_dir,ds$ids[kk]))
    V = VF[[1]]
    FF = VF[[2]]
	
    #Scale according to highest resolution point cloud
    V = V - repmat(ds$shape[[kk]]$center,1,dim(V)[2])
    if(scale){
      V = V / (ds$shape[[kk]]$scale / sqrt(ds$N[ds$K]) )
    }
    
    #Rotate according to ga
    if (det(ga$R[[kk]]) <0){ 
      FF = FF[c(2,1,3),]
      print("Determinant less than 0, flipping model.")
	  meshes.flipped<-c(meshes.flipped, as.character(ds$ids[kk]))
	  } 
    V = ga$R[[kk]]%*%V
    
	#Added additional outputs to get center and rotation matrix used for each bone
	write_center_path = file.path(output_dir,paste(substr(ds$ids[kk], 1, nchar(ds$ids[kk])-4),"_center.txt",sep=""))
	write.matrix(ds$shape[[kk]]$center, file = write_center_path, sep = " ")
	write_matrix_path = file.path(output_dir,paste(substr(ds$ids[kk], 1, nchar(ds$ids[kk])-4),"_rotation_matrix.txt",sep=""))
	write.matrix(ga$R[[kk]], file = write_matrix_path, sep = " ")
    write_off_path = file.path(output_dir,paste(substr(ds$ids[kk], 1, nchar(ds$ids[kk])-4),".off",sep=""))
    write_off(write_off_path,V,FF,FALSE)    
	
  }
  write.table(meshes.flipped, file = paste(output_dir, "/flipped_meshes.txt", sep=""), sep = ",")
}
 