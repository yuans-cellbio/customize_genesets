\name{read_gmt_term2name}
\alias{read_gmt_term2name}
\title{Read GMT Term Names}
\description{
Read the term identifier and description columns from a GMT file into a
\code{TERM2NAME}-style data frame.
}
\usage{
read_gmt_term2name(filepath)
}
\arguments{
  \item{filepath}{Path to a GMT file.}
}
\value{
A data frame with columns \code{ID} and \code{Name}. Empty GMT files return an
empty data frame with the same columns.
}
