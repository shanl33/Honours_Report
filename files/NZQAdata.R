nzqa2016 <- read.csv("http://www.nzqa.govt.nz/assets/Studying-in-NZ/Secondary-school-and-NCEA/stats-reports/2016/Qualification-Statistics-School-2016-29032017.csv")
# Drop vars that will not be used (eg. cumulative achievement)
nzqa2016 <- nzqa2016[,-c(1,8,10)]
names(nzqa2016) <- c("Decile", "Region", "School", "Year", "Qualification", "Achieve_participate", "Achieve_roll", "Small_sch")
levels(nzqa2016$Qualification) <- c("L1", "L2", "L3", "UE")
# Subset to use only Year 11 with Level 1, etc
nzqa <- nzqa2016 %>% filter(((Qualification=="L1") & (Year==11)) | 
                              ((Year==12) & (Qualification=="L2")) |
                              ((Year==13) & (Qualification=="L3")) | 
                              ((Year==13) & (Qualification=="UE"))) %>%
  # Reshape so that each row represents an observation (a school)
  # Current.Achievement.Rate.Participation kept for analysis only
  spread(Qualification, Achieve_participate, fill=0) %>%
  group_by(School) %>%
  summarise_at(c("L1", "L2", "L3", "UE"), sum) %>%
  inner_join(nzqa2016[, c(1, 2, 3, 8)]) %>% # Add Decile and Region and Small_sch variables
  distinct() %>% 
  filter(!((L1==0)&(L2==0)&(L3==0))) #Remove schools with 0% achievement rate for all levels

# Function to replace 0% with NA
zeros <- function(col) {
  replace(col, col==0, NA)
}

nzqa[2:5] <- sapply(nzqa[2:5], zeros)
# Remove schools with 'small cohort' warning (obscures pattern in non-small cohorts).
nzqa <- nzqa[nzqa$Small_sch=="",] # 439 school left
nzqa <- nzqa[, -8] # Remove Small_sch var
# Remove schools with any NA values 
nzqa <- nzqa[complete.cases(nzqa),] # 408 NZ schools
nzqa$Decile <- as.factor(nzqa$Decile)
# Subset Auckland schools only
akl <- as.data.frame(nzqa[nzqa$Region=="Auckland", -7]) # 91 AKL schools
rownames(akl) <- akl$School # Need for tooltips in tour
