#' Algorithm 1: SIGDROP candidate PSA declines
#'
#' Identify candidate PSA declines after ordering measurements by time.
#'
#' @param psa Numeric vector of nonmissing PSA values in ng/mL.
#' @param days Numeric vector of nonmissing measurement times in days, using a
#'   common origin (e.g., days since diagnosis). Must have the same length as
#'   `psa`.
#' @return A numeric matrix with one row per candidate decline and columns
#'   `drop_day` (peak time), `nadir_day` (first threshold-crossing time),
#'   `PSAmax` (peak value), and `PSAmin` (value at that threshold crossing).
#'   Returns `NULL` if no candidate is found or fewer than two values are supplied.
#' @examples
#' find_sigdrop(psa = c(8, 0.05), days = c(0, 30))
find_sigdrop <- function(psa, days) {
    stopifnot(length(psa) == length(days))
    # Ensure sorted by time (in case input isn't ordered)
    ord <- order(days)
    psa  <- psa[ord]
    days <- days[ord]

    n <- length(psa)
    
    if (n <= 1) return(NULL)
    
    # Initialize PSAmax
    PSAmax <- psa[1]
    day_PSAmax <- days[1]
    # PSAmin <- NA
    # drop_date <- NA
    # nadir_date <- NA
    
    drop_data <- c()

    for (j in 1:(n - 1)) {
            delta <- days[j + 1] - day_PSAmax
            psa_diff <- psa[j] - psa[j + 1]
            if (psa_diff <= 0) {
                # Update new maximum if
                # PSA increases or
                # more than 1 year has passed since last maximum
                if (psa[j + 1] >= PSAmax || delta > 365) {
                    PSAmax <- psa[j + 1]
                    day_PSAmax <- days[j + 1]
                }
            } else{
                beta <- PSAmax - psa[j + 1]
                alpha <- beta / max(PSAmax, 0.0001)
                if ((alpha >= 0.75 && beta >= 3) || (alpha >= 0.5 && beta >= 4)) {
                    PSAmin <- psa[j + 1]
                    drop_day <- day_PSAmax
                    nadir_day <- days[j + 1]
                    
                    item <- c(
                        drop_day = drop_day,
                        nadir_day = nadir_day,
                        PSAmax = PSAmax,
                        PSAmin = PSAmin)
                    drop_data <- rbind(drop_data, item)

                    # reset max after drop
                    # ignore consecutive drop
                    if (j+2 > n) break
                    day_PSAmax <- days[j + 2]
                    PSAmax <- psa[j + 2]
                }else{
                    if (delta> 365) {
                        PSAmax <- psa[j + 1]
                        day_PSAmax <- days[j + 1]
                    } else{
                        # look ahead guard, in case if the next 2 already passed one year and we have a big gap between max and next 2
                        # e.g.: max -- (300 days gap) -- next1 -- (300 days gap)- next2
                        if (j+2<= n){
                            gamma <- days[j+2] - day_PSAmax
                            if (gamma> 365) {
                                PSAmax <- psa[j + 1]
                                day_PSAmax <- days[j + 1]
                            }
                        }
                    }
                }
            }
    }
    if (length(drop_data)>0){
        return(drop_data)
    }else{
        return(NULL)
    }
}
