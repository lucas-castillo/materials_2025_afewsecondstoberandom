# Data
## survey_table.csv
Demographic questionnaire:
- **uuid**: participant id
- **age**: reported age in years
- **gender**: reported gender
- **genderdescr**: if participant chooses 'Other' as gender, they can describe
- **lifespan**: Answer to the question: which lifespan do you think you said the most

## rg_table.csv
- **uuid**: participant id
- **condition**: Fast (80 items per minute) or Slow (40 items per)
- **block**: Whether this is the 1st or 2nd second the participant produced
- **index**: Item index in the sequence
- **start**: start time of the vocalization, in seconds
- **end**: end time of the vocalization, in seconds
- **value**: lifespan produced
- **blockorder**: whether the participant did the Fast condition first (FS) or second (SF).
