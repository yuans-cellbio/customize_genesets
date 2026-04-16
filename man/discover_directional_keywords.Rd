\name{discover_directional_keywords}
\alias{discover_directional_keywords}
\title{Inspect Candidate Directional GO Keywords}
\description{
Scan GO Biological Process term names for directional wording and summarize
candidate keywords that may be useful when refining the signed GO keyword
lists.
}
\usage{
discover_directional_keywords(
  pos_keywords = DEFAULT_POS_KEYWORDS,
  neg_keywords = DEFAULT_NEG_KEYWORDS,
  candidate_words = NULL
)
}
\arguments{
  \item{pos_keywords}{Current positive-direction keyword list.}
  \item{neg_keywords}{Current negative-direction keyword list.}
  \item{candidate_words}{Optional character vector of candidate phrases to
  scan. Defaults to a broader built-in list of directional phrases.}
}
\value{
A data frame with one row per detected keyword and columns describing how many
GO BP terms matched, whether the keyword is already known, and example terms.
}
