GMS Statistics Programme wrap-up session
==============================


First we'll build a data frame of talks:
```R
talks = list(
	Arni = "Genomic inflation and genomic control",
	Jennifer = "Gibbs sampling",
	Lino = "Statistical methods for analysing hierarchical phenotype data",
	Minette = "Bayes Theorem",
	Samvida = "Mendelian randomisation",
	Nina = "Associating endometriosis and immunological diseases using population-scale health records",
	Emine = "Topic modelling and clustering for scATAC-seq data",
	Jiayuan = "Fine-mapping"
)

# Turn them into a data frame
talks = data.frame(
	who = names( talks ),
	title = unlist( talks ),
	type = c( rep( "GMS", 5 ), rep( "non-GMS", 3 )),
	row.names = NULL
)

```


Now we'll set the initial probabilities. In 2019 we have more GMS students than non-GMS students,
so let's pick a GMS student talk first:

```R
initialProbs = c( GMS = 1, 'non-GMS' = 0 )
```

We would like the HMM to prefer to alternate between GMS and non-GMS student talks. So let's give
the HMM a 7/8th chance of switching between types:

```
#                       prob = 7/8
#                      <----------
# GMS student talks                    non-GMS student talks
#                      ----------->
#                       prob = 7/8
```

Like this:
```R
transitions = matrix(
	nrow = 2, ncol = 2, dimnames = list(
		c( "GMS", "non-GMS" ),
		c( "GMS", "non-GMS" )
	)
)
transitions["GMS",] = c( 1/8, 7/8 )
transitions["non-GMS",] = c( 7/8, 1/8 )
```

Here is some code to run the HMM.  I like to try to collect related code into objects - I'll use an R list to collect the various functions.  (A better programmer might refer to the [R OO field guide](http://adv-r.had.co.nz/OO-essentials.html) here.)
```R
hmm = list(

	# sample(): a helper function to pick one item from a set of possibilities
	sample = function( from, probabilities ) {
		u = runif(1)
		return( from[min(which( cumsum(probabilities) >= u ))] )
	},

	# emit(): pick a talk uniformly among remaining talks with the given type
	# if there's no such talk left we return an 'empty' talk (i.e. NULL)
	emit = function( talks, schedule, type ) {
		w = which( talks$type == type & !talks$who %in% schedule$who )
		if( length(w) > 0 ) {
			pick = hmm$sample( w, probabilities = rep( 1/length(w), length(w) ))
			return( talks[pick, ] )
		} else {
			return(NULL)
		}
	},

	# initialise(): Pick an initial talk
	initialise = function( talks, initialProbs ) {
		schedule = data.frame()
		talkType = hmm$sample( names( initialProbs ), probabilities = initialProbs )
		talk = hmm$emit( talks, schedule, talkType )
		schedule = rbind( schedule, talk )
		return( schedule )
	},

	# step(): Pick the next talk, given last currently scheduled talk type.
	step = function( talks, schedule, transitions ) {
		currentTalk = schedule[nrow( schedule ), ]
		# transition to next type of talk
		nextTalkType = hmm$sample( c( "GMS", "non-GMS" ), probabilities = transitions[currentTalk$type,] )
		# pick a new talk
		nextTalk = hmm$emit( talks, schedule, nextTalkType )
		if( is.null( nextTalk )) {
			cat( sprintf( "No %s talks left to emit, emitting an empty talk.\n", nextTalkType ))
		} else {
			schedule = rbind( schedule, nextTalk )
			cat( sprintf( "Emitted the talk from %s.  Schedule is:\n", nextTalk$name ))
		}
		return( schedule )
	}
)
```

We're set! Let's run it:
```R
	schedule = hmm$initialise( talks, initialProbs ); print( schedule ) ;
	while( nrow( schedule ) < nrow( talks )) {
		schedule = hmm$step( talks, schedule, transitions ); print( schedule ) ;
	}
	print( schedule )
```

Finally, let's add times and so on: let's print the full schedule:
```R
formatSchedule <- function(
	schedule,
	startTimeInMinutesPastMidnight = (12*60) + 20,
	talkDurationInMinutes = 5
) {
	# Function to print a talk start and end time
	formatTime <- function( start.time, duration = 5 ) {
		end.time = start.time + duration
		start.hour = floor( start.time / 60 )
		end.hour = floor( end.time / 60 )

		start.minute = start.time - 60 * start.hour
		end.minute = end.time - 60 * end.hour

		sprintf( "%.2d:%.2d - %.2d:%.2d", start.hour, start.minute, end.hour, end.minute )
	}
	times = seq( from = startTimeInMinutesPastMidnight, length = nrow( schedule ), by = talkDurationInMinutes + 2 )
	schedule$time = formatTime( times )
	schedule = rbind(
		data.frame( who = "-", title = "Lunch / setup", type = "-", time = "12:00 - 12:20" ),
		schedule,
		data.frame( who = "-", title = "Wrap-up", type = "-", time = "13:15 - 13:30" )
	)
	return( schedule )
}

print( formatSchedule( schedule ))
```

