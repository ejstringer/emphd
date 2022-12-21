
# setup ------------------------------------------------------------------------

#' setup repository
#' 
#' Setup the folders and git ignore files at the beginning of a repo.
#' Folder structure is so
#' 
#' * code (pushed)
#' * data (ignored)
#' * figures (ignored)
#' * output (ignored)
#' 
#' @param pushFolder - folders that will be pushed to github
#' @param ignoreFolder - folders whose contents will not be pushed to github
#' @return creates folders and .gitignore files
#' @export
#' @author Emily Stringer
#' @examples em.setup.repo()
#' 
#' ### make function in emphd !! ###
#' 

em.setup.repo <- function(pushFolders = "code", 
                          ignoreFolder = c("data", "figures", "output")){
folders <- c(pushFolders, ignoreFolder)

lapply(folders, dir.create)

ignore <- c("# Ignore everything in this directory",
            "*", 
            "# Except this file",
            "!.gitignore")

ignoreFolderPath <- paste0("./", ignoreFolder, "/.gitignore")

lapply(ignoreFolderPath, function(x) write(ignore, x))



bigIgnore <- ".Rhistory
              .Rapp.history
                     .RData
                 .Ruserdata
                     *-Ex.R
                  /*.tar.gz
                 /*.Rcheck/
               .Rproj.user/
           vignettes/*.html
          vignettes/*.pdf
              .httr-oauth
                  *_cache/
                  /cache/
                 *.utf8.md
                 *.knit.md
                 .Renviron
                    *.html
                    *.docx
                     *.pdf
                     *.png
                     *.jpg"
bigIgnore2 <- gsub(" ", "", bigIgnore)
bigIgnore3 <- gsub("[\r\n]", " ", bigIgnore2)
bigIgnore4 <- do.call("c", strsplit(bigIgnore3, " "))

write(bigIgnore4, file = "./.gitignore")

}