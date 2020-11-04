# Utility Functions
# Author: Koki Sasagawa 

create_file_directory = function(file_path){
  if(!dir.exists(file_path)) {
    dir.create(file_path)
  }
}