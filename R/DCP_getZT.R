#' Title Calculate Zeitgeber Time
#'
#' Calculate Zeitgeber Time based on clock time based on local sunrise time.
#'
#' @param t Recorded time in POSIXlt format. Please make sure that the input time has the correct time zone.
#' @param lat Local latitude corresponding to the recorded time.
#' @param long Local longitude corresponding to the recorded time.
#' @param ZT.min The mininum of output ZT. The default is -6, which corresponds to output range of (-6, 18).
#'
#' @return Zeitgeiber time corrected by local sunrise
#' @export
#'
#' @examples
#' t = as.POSIXlt("2022-09-17 13:43:00 EST", tz="America/New_York",usetz=TRUE)
#' lat = 40.4406; long = -79.9959 #This is the Pittsburgh coordinates
#' t.correct = DCP_getZT(t, lat, long, ZT.min = -6)
DCP_getZT = function(t, lat, long, ZT.min = -6){

  message("Please make sure that the input time has the correct time zone. ")
  #check if input is POSIXlt
  stopifnot("Please input t as POSIXlt format. \n
            Time data can be converted to POSIXlt with as.POSIXlt(). Notice that Sys.timezone() will be used if you did not specify correct tz in as.POSIXlt(). " =
              "POSIXlt" %in% class(t))
  stopifnot("lat and long should both be vectors" = is.vector(lat)&is.vector(long))
  stopifnot("Please make sure that t, lat, and long are of the same lengths" = length(t)==length(lat)&length(t)==length(long))

  #calculate sunrise time
  day = t$yday
  prev_sunrise<- SunRiseTime(day-1, Lat = lat, Long = long)
  cur_sunrise <- SunRiseTime(day, Lat = lat, Long = long)
  next_sunrise <- SunRiseTime(day+1, Lat = lat, Long = long)

  t_chr = as.character(t)
  t_UTC = lubridate::with_tz(t, tzone='UTC') #convert time to UTZ
  t_UTC_num = t_UTC$hour + t_UTC$min/60
  # t_num = t$hour+t$min/60
  ZT.max = ZT.min+24


  date.time.flag = 0
  corrected.time = vector()
  for(i in seq_along(t)){
    a.t = t[i]
    if(t_UTC_num[i]+cur_sunrise$timezone[i]-cur_sunrise$sunrise[i]>=0){
      corrected.time[i]<- t_UTC_num[i]+cur_sunrise$timezone[i]-cur_sunrise$sunrise[i]
    } else {
      corrected.time[i] <- t_UTC_num[i]+prev_sunrise$timezone[i] + (24 - prev_sunrise$sunrise[i])
    }

    #transform the time into (ZT.min, ZT.max)
    ifelse(corrected.time[i]>ZT.max, corrected.time[i]-24,
           ifelse(corrected.time[i]<ZT.min, corrected.time[i]+24, corrected.time[i]))

    #check if t contains both date and time
    if(is.na(strptime(a.t, format="%Y-%m-%d %H:%M:%S"))){
      date.time.flag = 1
      corrected.time[i] = NA
    }
  }
  if(date.time.flag){
    warning("t without time (HH:MM:SS) input will be set to NA. ")
  }

  return(corrected.time)
}

SunRiseTime = function(d,Lat,Long){
  ## d is the day of year
  ## Lat is latitude in decimal degrees
  ## Long is longitude in decimal degrees (negative == West)

  ##This method is copied from:
  ##Teets, D.A. 2003. Predicting sunrise and sunset times.
  ##  The College Mathematics Journal 34(4):317-321.

  ## At the default location the estimates of sunrise and sunset are within
  ## seven minutes of the correct times (http://aa.usno.navy.mil/data/docs/RS_OneYear.php)
  ## with a mean of 2.4 minutes error.

  ## Function to convert degrees to radians
  rad<-function(x)pi*x/180

  ##Radius of the earth (km)
  R=6378

  ##Radians between the xy-plane and the ecliptic plane
  epsilon=rad(23.45)

  ##Convert observer's latitude to radians
  L=rad(Lat)

  ## Calculate offset of sunrise based on longitude (min)
  ## If Long is negative, then the mod represents degrees West of
  ## a standard time meridian, so timing of sunrise and sunset should
  ## be made later.
  timezone = -4*(abs(Long)%%15)*sign(Long)
  timezone2 = abs(Long)%/%15*sign(Long)

  ## The earth's mean distance from the sun (km)
  r = 149598000

  theta = 2*pi/365.25*(d-80)

  z.s = r*sin(theta)*sin(epsilon)
  r.p = sqrt(r^2-z.s^2)

  t0 = 1440/(2*pi)*acos((R-z.s*sin(L))/(r.p*cos(L)))

  ##a kludge adjustment for the radius of the sun
  that = t0+5

  ## Adjust "noon" for the fact that the earth's orbit is not circular:
  n = 720-10*sin(4*pi*(d-80)/365.25)+8*sin(2*pi*d/365.25)

  ## now sunrise and sunset are:
  sunrise = (n-that+timezone)/60
  sunset = (n+that+timezone)/60

  return(list("sunrise" = sunrise,"sunset" = sunset, "timezone" = timezone2))
}
