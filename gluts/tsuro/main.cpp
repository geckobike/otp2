/*

	Pre place cards, extra points for using cards
	Line length for scoring
	Rethink Undo and card dropping
	Player position
	Slow moving timer mode - cards that are slow and faster

*/

#include <time.h>
#include <sys/time.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "engine/engine.h"
#include <assert.h>

#include <GL/gl.h>
#include <GL/glut.h>

#ifdef _MSC_VER
#pragma warning(disable:4244 4305)  // for VC++, no precision loss complaints
#endif

#define clampf(x,min,max) ((x) < (min)) ? (min) : ( (x) > (max) ? (max) : x )

static int s_width = 480 * 2;
static int s_height = 320 * 2;

static bool s_paused = false;
static float s_averageDt = 0.f;

#define XYZ(v) (v).x, (v).y, (v).z
#define XYZp(v) (v)->x, (v)->y, (v)->z

static void drawCircle(float x, float y, float r, bool filled = true);

//======================================================================================================================================================
//======================================================================================================================================================

//======================================================================================================================================================
//======================================================================================================================================================

struct Timer
{
	timeval time;
};

void timerUpdate(Timer* timer)
{
	gettimeofday(&timer->time,NULL);
}

double timerGetTickSinceLastUUpdate(Timer* timer)
{
    struct timeval now;
 
    long long sec;
    long long microsec;
    
    gettimeofday(&now,NULL);
    
    sec = now.tv_sec - timer->time.tv_sec;
    microsec = now.tv_usec - timer->time.tv_usec;
   
    return (double)sec + ( (double)microsec / (1e6) );
}

static Timer g_time;
static Timer g_demoTimer;

static void dlDrawLine(const vec3* from, const vec3* to)
{
	glBegin(GL_LINES);
	glVertex3f(from->x, from->y, from->z);
	glVertex3f(to->x, to->y, to->z);
	glEnd();
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////

/*
///////////////////////////////////////////////////////////////////////////////////////////////////////////

	Card
          5  4
       ___|__|___
      |          |
   6 -|          |- 3
      |          |
   7 -|          |- 2
      |__________|
          |  |
          0  1

       ___________
      |    |  |   |
      |_    --   _|   
      | |       | |  CARD 1
      |_|       |_|
      |     __    |
      |____|__|___|
    
       ___________
      |    |  |   |
      |___ |   \__|   
      |    \\     |  CARD 2
      |___  \\____|
      |   \  \    |
      |____|__|___|
    
       ___________
      |    |  |   |
      |___  \_|___|   
      |   \   |   |  CARD 3
      |____|_  \__|
      |    | \    |
      |____|__|___|
    
       ___________
      |    |  |   |
      |_ /     \ _|   
      |           |  CARD 4
      |_         _|
      |  \     /  |
      |____|__|___|
    
       ___________
      |    |  |   |
      |____|__|___|   
      |    |  |   |  CARD 5
      |____|__|___|
      |    |  |   |
      |____|__|___|
       ___________
      |    |  |   |
      |____|  |___|   
      |    \\//   |  CARD 6
      |____//\\___|
      |    |  |   |
      |____|__|___|

       ___________
      |    |  |   |
      |____|__|___|   
      |    \__/   |  CARD 7
      |____/__\___|
      |    |  |   |
      |____|__|___|


       ___________
      |    |  |   |
      |__ /    \__|   
      |           | CARD 8
      |_____  ____|
      |     /\    |
      |____|__|___|

       ___________
      |    |  |   |
      |__ / __|_ _|   
      |    /  |   | CARD 9
      |____|_  \__|
      |    | \    |
      |____|__|___|

       ___________
      |    |  |   |
      |____|  |_ _|   
      |    /\/|   | CARD 10 (drawing not so good!)
      |___/ /\ \__|
      |    |  \   |
      |____|__|___|

       ___________
      |    |__|   |
      |___________|   
      |           | CARD 11
      |___________|
      |     __    |
      |____|__|___|

       ___________
      |    |  |   |
      |_         _|   
      |           |
      |_         _|
      |           |
      |____|__|___|

       ___________
      |    |  |   |
      |_         _|   
      |           |
      |_         _|
      |           |
      |____|__|___|

       ___________
      |    |  |   |
      |_         _|   
      |           |
      |_         _|
      |           |
      |____|__|___|


///////////////////////////////////////////////////////////////////////////////////////////////////////////
*/

#define MOUSEUP 1
#define MOUSEDOWN 2
#define MOUSEDRAG 4

const int tsuroNumCells = 6;
const int tsuroVectorPathHalfSize = 12;
const int tsuroVectorPathSize = (tsuroVectorPathHalfSize*2 + 1);
const float tsuroCardSize = 96.f;
const float tsuroCardSizePad = 1.f;
int tsuroNumBaseCards = 35;
float tsuroPlayerSpeed = 4.0f;

// Helpers for determining which edges the paths connect
static const int tsuroEdgeOpposite1[8] = {5,4,7,6,1,0,3,2};	// Opposite edge (directly)
static const int tsuroEdgeOpposite2[8] = {4,5,6,7,0,1,2,3};	// Opposite edge
static const int tsuroEdgeSame[8] = {1,0,3,2,5,4,7,6};	// Same edge

/*
	The tile transition for each path that can be taken
*/
static const int tsuroTileTransition[8][2]=
{
	{0,-1},		// If we end at path zero, then move down one tile
	{0,-1},

	{+1,0},
	{+1,0},

	{0,+1},
	{0,+1},

	{-1,0},
	{-1,0},
};


static const float tsuroCardPoints[8][2] = 
{
	{-tsuroCardSize*(1.f/6.f), -0.5f*tsuroCardSize},
	{+tsuroCardSize*(1.f/6.f), -0.5f*tsuroCardSize},

	{+0.5f*tsuroCardSize, -tsuroCardSize*(1.f/6.f)},
	{+0.5f*tsuroCardSize, +tsuroCardSize*(1.f/6.f)},

	{+tsuroCardSize*(1.f/6.f), +0.5f*tsuroCardSize},
	{-tsuroCardSize*(1.f/6.f), +0.5f*tsuroCardSize},

	{-0.5f*tsuroCardSize, +tsuroCardSize*(1.f/6.f)},
	{-0.5f*tsuroCardSize, -tsuroCardSize*(1.f/6.f)},
};

//===================================================
// Vector path used for drawing and route animation
//===================================================
struct tsuroVectorPath
{
	vec3 centre[tsuroVectorPathSize];	// central path - used for the animation
	vec3 inner[tsuroVectorPathSize];		// inner part of track
	vec3 outer[tsuroVectorPathSize];		// outer part of track
};

struct tsuroCard
{
	int paths[8];
	tsuroVectorPath vpaths[8];
	vec3 colour[8];
};

struct tsuroNode
{
	tsuroCard card;
	int used;
};

struct tsuroPlayer
{
	struct State
	{
		int x,y,path;		// tile x,y and path
	};

	State start;			// State at the start of the game
	State current;			// Current State
	State target;			// State the player will reach after a card has been placed.
		
	float animFraction;		// Current fraction through the tile
	float totalLength;		// Total distance travelled

	vec3 colour;
};

struct tsuroGame
{
	tsuroNode nodes[tsuroNumCells][tsuroNumCells];
	tsuroCard deck[tsuroNumCells*tsuroNumCells-1];
	tsuroCard selection[3];
	tsuroPlayer player;
	int cardIndex;				// Index into original deck
	bool selectionUsed[3];
	int mode;
	int state;

	// Undo
	struct State
	{
		tsuroPlayer player;
		tsuroCard selection[3];
		bool selectionUsed[3];
		int cardIndex;
	};
	State last;
};


#define TSURO_MODE_SINGLE_PLAYER	0
#define TSURO_MODE_DEMO				1
#define TSURO_MODE_NUM				2

#define TSURO_STATE_GAME_OVER		0
#define TSURO_STATE_PLAYING			1

///////////////////////////////////////

static void tsuroMovePlayer(float frac);
static void tsuroMouseButton(int state, int x, int y);

///////////////////////////////////////

static void generateQuarterCurve(vec3* out, const vec3* p, const vec3* x, const vec3* y, int N, float gamma=1.f)
{
	if (N>1)
	{
		const float dAlpha = PI * 0.5f / float(N-1);
		float alpha = 0.f;
		vec3* v = out;
		vecadd(v, p, x);
		for (int i=1; i<N; i++)
		{
			v++;
			alpha += dAlpha;
			vecaddscale(v, p, x, powf(cosf(alpha),gamma));
			vecaddscale(v, v, y, powf(sinf(alpha),gamma));
		}
	}
}

void tsuroCardDump(tsuroCard* card)
{
	printf("dump card!\n");
	for (int i=0; i<8; i++)
	{
		printf("%d -> %d, \n", i, card->paths[i]);
	}
	printf("\n");
}

static void tsuroCardSetAllColours(tsuroCard* card, float r, float g, float b)
{
	for (int i=0; i<8; i++)
	{
		vecset(&card->colour[i], r, g, b);
	}
}

void tsuroCardCreate(tsuroCard* card, int input[4][2])
{
	bool okay = true;
	memset(card, 0, sizeof(card));
	memset(card->paths, 0xff, sizeof(card->paths));

	tsuroCardSetAllColours(card, 1.f, 1.f, 1.f);

	for (int i=0; i<4; i++)
	{
		int from = input[i][0];
		int to = input[i][1];
		if ((card->paths[from] & card->paths[to]) == -1)
		{
			card->paths[from] = to;
			card->paths[to] = from;
		}
		else
		{
			printf("WARNING: Invalid input for card!\n");
			okay = false;
		}
	}

	if (okay)
	{
		// Generate the vector path
		const float s = 0.5f*tsuroCardSize;				// half width/scale of card
		const float a = tsuroCardSize*(1.f/6.f);		// position of connection points

		static const float points[8][2] = 
		{
			{-a, -s},
			{+a, -s},

			{+s, -a},
			{+s, +a},

			{+a, +s},
			{-a, +s},

			{-s, +a},
			{-s, -a},
		};

		for (int i=0; i<8; i++)	// We are doing this twice (but its easier this way!)
		{
			int from = i;
			int to = card->paths[i];
			vec3 f = {points[from][0], points[from][1], 0.f};	// from
			vec3 t = {points[to][0], points[to][1], 0.f};		// to 

			if (to == tsuroEdgeOpposite1[from])
			{
				// Directly opposite
				vec3 dv;
				vecsub(&dv, &t, &f);
				vecscale(&dv, &dv, 1.f/((float)(tsuroVectorPathSize-1)));
				vec3* v = card->vpaths[i].centre;
				v[0] = f;
				for (int n=1; n<tsuroVectorPathSize; n++)
				{
					vecadd(&v[n], &v[n-1], &dv);
				}
			}
			else if (to == tsuroEdgeOpposite2[from])
			{
				// Opposite wall

				vec3 centre = {0.f, 0.f, 0.f};

				vec3 adjacentFrom = {points[tsuroEdgeSame[from]][0], points[tsuroEdgeSame[from]][1], 0.f};	// The point that is on the same edge as from
				vec3 adjacentTo = {points[tsuroEdgeSame[to]][0], points[tsuroEdgeSame[to]][1], 0.f};			// The point that is on the same edge as to

				vec3 focal1, focal2;
				vecmidpoint(&focal1, &adjacentFrom, &f);
				vecmidpoint(&focal2, &adjacentTo, &t);

				vec3 x,y;

				vecsub(&x, &centre, &focal1);
				vecsub(&y, &f, &focal1);
				generateQuarterCurve(&card->vpaths[i].centre[0], &focal1,&y,&x, tsuroVectorPathHalfSize+1, 1.5f);

				vecsub(&x, &centre, &focal2);
				vecsub(&y, &t, &focal2);
				generateQuarterCurve(&card->vpaths[i].centre[tsuroVectorPathHalfSize], &focal2,&x,&y, tsuroVectorPathHalfSize+1, 1.5f);
			}
			else if (to == tsuroEdgeSame[from])
			{
				vec3 centre = {0.f, 0.f, 0.f};

				// Same wall
				vec3 midpoint;
				vecmidpoint(&midpoint, &f, &t);

				vec3 x,y;
				vecsub(&x, &f, &midpoint);
				vecsub(&y, &centre, &midpoint);
				vecscale(&y, &y, 0.5f);
				generateQuarterCurve(&card->vpaths[i].centre[0], &midpoint,&x,&y, tsuroVectorPathHalfSize+1, 0.8);
				vecneg(&x, &x);
				generateQuarterCurve(&card->vpaths[i].centre[tsuroVectorPathHalfSize], &midpoint,&y,&x, tsuroVectorPathHalfSize+1, 0.8);
			}
			else
			{
				vec3 focal;
				switch(from)
				{
					case 0: focal.y = -s; break;
					case 1: focal.y = -s; break;
					case 2: focal.x = +s; break;
					case 3: focal.x = +s; break;
					case 4: focal.y = +s; break;
					case 5: focal.y = +s; break;
					case 6: focal.x = -s; break;
					case 7: focal.x = -s; break;
				}
				switch(to)
				{
					case 0: focal.y = -s; break;
					case 1: focal.y = -s; break;
					case 2: focal.x = +s; break;
					case 3: focal.x = +s; break;
					case 4: focal.y = +s; break;
					case 5: focal.y = +s; break;
					case 6: focal.x = -s; break;
					case 7: focal.x = -s; break;
				}

				vec3 x,y;
				vecsub(&x, &f, &focal);
				vecsub(&y, &t, &focal);
				generateQuarterCurve(&card->vpaths[i].centre[0], &focal,&x,&y, tsuroVectorPathSize, 1.2f);

				// printf("DUMP: from = %f %f %f\n", XYZ(f));
				// vec3* v = card->vpaths[i].centre;
				// for (int n=0; n<tsuroVectorPathSize; n++)
				// {
				// 	printf("%f %f %f\n", XYZp(v));
				// 	v++;
				// }
			}
		}
	}
}

static void tsuroCardCopy(tsuroCard* dst, tsuroCard* src)
{
	memcpy(dst, src, sizeof(*src));
}

static void tsuroCardSwap(tsuroCard* a, tsuroCard* b)
{
	tsuroCard tmp;
	tsuroCardCopy(&tmp, a);
	tsuroCardCopy(a, b);
	tsuroCardCopy(b, &tmp);
}

static void rotateVectors(vec3* v, int count, int rotate)
{
	switch(rotate)
	{
		case 1:
			for (int i=0; i<tsuroVectorPathSize; i++)
			{
				float x = v[i].x;
				float y = v[i].y;
				v[i].x = -y;
				v[i].y = x;
			}
			break;
		case 2:
			for (int i=0; i<tsuroVectorPathSize; i++)
			{
				float x = v[i].x;
				float y = v[i].y;
				v[i].x = -x;
				v[i].y = -y;
			}
			break;
		case 3:
			for (int i=0; i<tsuroVectorPathSize; i++)
			{
				float x = v[i].x;
				float y = v[i].y;
				v[i].x = y;
				v[i].y = -x;
			}
			break;
	}
}

void tsuroCardRotate(tsuroCard* card, int rotate)
{
	rotate = (rotate + 4) % 4;

	// Rotate the drawing
	for (int p=0; p<8; p++)
	{
		rotateVectors(card->vpaths[p].centre, tsuroVectorPathSize, rotate);
	}

	rotate = (rotate * 2) % 8;
	
	tsuroVectorPath tmpPath[8];
	for (int p=0; p<8; p++)
	{
		memcpy(tmpPath[p].centre, card->vpaths[p].centre, sizeof(card->vpaths[p].centre));
	}
	for (int p=0; p<8; p++)
	{
		memcpy(card->vpaths[(p+rotate)%8].centre, tmpPath[p].centre, sizeof(card->vpaths[p].centre));
	}

	int tmp[8];
	memcpy(tmp, card->paths, sizeof(card->paths));

	for (int i=0; i<8; i++)
	{
		card->paths[(i+rotate)%8] = (tmp[i] + rotate)%8;
	}
}

static inline float tsuroCardScreenX(int i)
{
	// Could probably do this nicer!
	float dx = tsuroCardSize + tsuroCardSizePad;
	float width= dx * tsuroNumCells;
	float x0 = tsuroCardSize*0.75f;
	return x0 + (float)i*dx;
}

static inline float tsuroCardScreenY(int j)
{
	return tsuroCardScreenX(j);
}

void tsuroCardDraw(tsuroCard* card, float x, float y)
{
	const float s = 0.5f*tsuroCardSize;
	const float eps = s * 0.1f;
	const float eps2 = s * 0.04f;

	// Draw border and main card
	glDisable(GL_POINT_SMOOTH);

	glColor3f(0.4f,0.3f,0.3f);
	glBegin(GL_QUADS);
		glVertex3f(x+eps2-s,y+eps2-s,0.f);
		glVertex3f(x-eps2+s,y+eps2-s,0.f);
		glVertex3f(x-eps2+s,y-eps2+s,0.f);
		glVertex3f(x+eps2-s,y-eps2+s,0.f);
	glEnd();

	glColor3f(1.f,0.f,0.f);
	glBegin(GL_LINE_STRIP);
		glVertex3f(x+eps-s,y+eps-s,0.f);
		glVertex3f(x-eps+s,y+eps-s,0.f);
		glVertex3f(x-eps+s,y-eps+s,0.f);
		glVertex3f(x+eps-s,y-eps+s,0.f);
		glVertex3f(x+eps-s,y+eps-s,0.f);
	glEnd();

			
	// Draw Paths
	int drawn[8] = {0};

	for (int p=0; p<8; p++)
	{
		int from = p;
		int to = card->paths[p];
		if (drawn[from]==0)
		{
			const vec3* c = &card->colour[p];
			glColor3f(c->x, c->y, c->z);
			drawn[from]=1;
			drawn[to]=1;
			glBegin(GL_LINE_STRIP);
			for (int n=0; n<tsuroVectorPathSize; n++)
			{
				const vec3* v = &card->vpaths[p].centre[n];
				glVertex3f(x+v->x, y+v->y, 0.f);
			}
			glEnd();
		}
	}
}

static bool g_paused = false;
static int g_rotate = 0;
static int g_use = 0;			// Begin using the card
static int g_used = 0;			// The card has been used
static int g_using = 0;			// The card is being used
static int g_select = 0;
static int g_mouseX = 0;
static int g_mouseY = 0;
static bool g_dragging = false;
static tsuroGame g_game;
	
static void tsuroShuffleDeck()
{
	for (int repeat=0; repeat<64; repeat++)
	{
		int i = rand() % tsuroNumBaseCards;
		int j = rand() % tsuroNumBaseCards;
		tsuroCardSwap(&g_game.deck[i], &g_game.deck[j]);
	}
}

static void tsuroDealCard(int i)
{
	if (g_game.cardIndex < tsuroNumBaseCards)
	{
		g_game.selection[i] = g_game.deck[g_game.cardIndex];
		g_game.selectionUsed[i] = true;
		g_game.cardIndex++;
	}
	else
	{
		g_game.selectionUsed[i] = false;
	}
}

static void tsuroSnapShot(tsuroGame::State* state)
{
	state->player = g_game.player;	// Keep note of the target position
	state->selection[0] = g_game.selection[0];
	state->selection[1] = g_game.selection[1];
	state->selection[2] = g_game.selection[2];
	state->selectionUsed[0] = g_game.selectionUsed[0];
	state->selectionUsed[1] = g_game.selectionUsed[1];
	state->selectionUsed[2] = g_game.selectionUsed[2];
	state->cardIndex = g_game.cardIndex;
}

static void tsuroRestore()	// Undo
{
	tsuroGame::State* last = &g_game.last;

	g_game.selection[0] = last->selection[0];
	g_game.selection[1] = last->selection[1];
	g_game.selection[2] = last->selection[2];
	g_game.selectionUsed[0] = last->selectionUsed[0];
	g_game.selectionUsed[1] = last->selectionUsed[1];
	g_game.selectionUsed[2] = last->selectionUsed[2];
	g_game.cardIndex = last->cardIndex;

	// Remove the last placed card (which is situated at the last target)
	g_game.nodes[last->player.target.x][last->player.target.y].used = false;

	// Reset the colours on the board and move the player from the start
	// to the last position
	for (int i=0; i<tsuroNumCells; i++)
	{
		for (int j=0; j<tsuroNumCells; j++)
		{
			tsuroNode* node = &g_game.nodes[i][j];
			tsuroCardSetAllColours(&node->card, 1.f, 1.f, 1.f);
		}
	}

	tsuroPlayer* player = &g_game.player;

	// Store the current location
	tsuroPlayer::State current = player->current;
	float animFraction = player->animFraction;
	float totalLength = player->totalLength;

	// Recolour the board and move the player as far as the current location
	player->totalLength = 0.f;
	player->animFraction = 0.f;
	player->current.x = player->start.x;
	player->current.y = player->start.y;
	player->current.path = player->start.path;

	while (1)
	{
		tsuroNode* node = &g_game.nodes[player->current.x][player->current.y];
		if (player->totalLength >= totalLength)
		{
			printf("breaking cause player->totalLength = %f, but totalLength = %f\n", player->totalLength, totalLength);
			player->totalLength = totalLength;
			player->animFraction = animFraction;
			player->current = current;
			break;
		}
		if (!node->used)
		{
			printf("breaking cause %d %d is empty\n", player->current.x, player->current.y);
			break;
		}
		tsuroMovePlayer(1.f);
	}

	player->target = last->player.target;

	g_game.state = TSURO_STATE_PLAYING;
}

static void tsuroReset()
{
	tsuroPlayer* player = &g_game.player;

	// Set player colour
	vecset(&player->colour, 1.f, 1.f, 0.f);

	switch(g_game.mode)
	{
		case TSURO_MODE_DEMO:
		{
			g_game.cardIndex = 0;
			player->current.x = 0;
			player->current.y = 3;
			player->current.path = 7;

			tsuroShuffleDeck();

			// Place
			int n=0;
			for (int j=0; j<tsuroNumCells; j++)
			{
				for (int i=0; i<tsuroNumCells; i++)
				{
					if (n<tsuroNumBaseCards)
					{
						tsuroNode* node = &g_game.nodes[i][j];
						node->used = 1;
						node->card = g_game.deck[n];
						n++;
					}
				}
			}
			break;
		}

		case TSURO_MODE_SINGLE_PLAYER:
		{
			player->current.x = 0;
			player->current.y = 3;
			player->current.path = 7;
	
			g_game.cardIndex = 0;
			for (int i=0; i<3; i++)
			{
				tsuroDealCard(i);
			}

			for (int j=0; j<tsuroNumCells; j++)
			{
				for (int i=0; i<tsuroNumCells; i++)
				{
					g_game.nodes[i][j].used = 0;
				}
			}
			break;
		}
	}

	player->animFraction = 0.f;
	player->totalLength = 0.f;
	player->start = player->current;
	player->target = player->current;

	g_game.state = TSURO_STATE_PLAYING;
	tsuroSnapShot(&g_game.last);

	tsuroShuffleDeck(); 

	tsuroMouseButton(MOUSEUP,1000,1000);
	g_use = 0;
	g_used = 0;
	g_using = 0;
	g_select = 1;
	g_dragging = 0;
	g_rotate = 0;
	g_select = -1;
}

static void tsuroChangeMode()
{
	g_game.mode = (g_game.mode + 1) % TSURO_MODE_NUM;
	tsuroReset();
}

static void tsuroCreateDeck()
{
	tsuroCard cards[tsuroNumBaseCards];
	int n;

#define CARD(a,b, c,d, e,f, g,h) { int _data[4][2] = {{a,b},{c,d},{e,f},{g,h}}; tsuroCardCreate((&cards[n]),_data); n++; } 

	n = 0;
	CARD( 0,1, 2,3, 4,5, 6,7);
	CARD( 0,7, 1,6, 2,5, 3,4);
	CARD( 0,6, 1,7, 2,4, 3,5);
	CARD( 0,7, 1,2, 3,4, 5,6);
	CARD( 0,5, 1,4, 2,7, 3,6);

	CARD( 0,4, 1,5, 2,6, 3,7);	// Need better drawing
	CARD( 0,4, 1,5, 2,7, 3,6);
	CARD( 1,7, 0,2, 3,4, 5,6);
	CARD( 1,7, 0,3, 2,4, 5,6);
	CARD( 1,6, 0,3, 2,4, 5,7);

	CARD( 2,7, 3,6, 0,1, 4,5);
	CARD( 3,7, 2,6, 0,1, 4,5);
	CARD( 3,6, 4,5, 0,2, 1,7);
	CARD( 3,4, 5,6, 0,2, 1,7);
	CARD( 3,5, 4,6, 0,2, 1,7);

	CARD( 2,5, 4,7, 0,3, 1,6);
	CARD( 0,7, 1,6, 2,4, 3,5);
	CARD( 2,7, 0,6, 1,3, 4,5);
	CARD( 0,3, 7,4, 1,5, 6,2);	// Need better drawing
	CARD( 2,7, 1,6, 0,3, 4,5);

	CARD( 2,7, 1,3, 0,4, 6,5);
	CARD( 2,7, 0,6, 1,5, 3,4);
	CARD( 1,5, 7,3, 0,2, 4,6);
	CARD( 0,4, 6,2, 1,7, 5,3);
	CARD( 0,3, 1,2, 4,5, 6,7);

	CARD( 0,3, 1,5, 4,2, 6,7);
	CARD( 5,2, 4,0, 1,3, 6,7);
	CARD( 5,3, 4,0, 1,2, 6,7);
	CARD( 1,2, 0,6, 7,4, 5,3);
	CARD( 0,6, 1,7, 2,3, 4,5);

	// Replicate the last 5 cards!
	CARD( 0,3, 1,5, 4,2, 6,7);
	CARD( 5,2, 4,0, 1,3, 6,7);
	CARD( 5,3, 4,0, 1,2, 6,7);
	CARD( 1,2, 0,6, 7,4, 5,3);
	CARD( 0,6, 1,7, 2,3, 4,5);

	if (n>tsuroNumBaseCards)
	{
		printf("error: you need to increase numBaseCards!\n");
		exit(0);
	}

	// Copy the deck
	memcpy(g_game.deck, cards, sizeof(g_game.deck));
}

static void tsuroStart()
{
	timerUpdate(&g_time);
	timerUpdate(&g_demoTimer);

	memset(&g_game, 0, sizeof(g_game));

	//g_game.mode = TSURO_MODE_DEMO;
	g_game.mode = TSURO_MODE_SINGLE_PLAYER;

	tsuroCreateDeck();
	tsuroReset();
}

static void tsuroDemoUpdate(float dt)
{
	tsuroPlayer* player = &g_game.player;

	if (!g_paused)
	{
		tsuroMovePlayer(dt * tsuroPlayerSpeed);
		float timeRunning = timerGetTickSinceLastUUpdate(&g_demoTimer);
		if (timeRunning>5.f)
		{
			tsuroReset();
			timerUpdate(&g_demoTimer);
		}
	}
	if (g_rotate)
	{

		for (int i=0; i<tsuroNumCells; i++)
		{
			for (int j=0; j<tsuroNumCells; j++)
			{
				tsuroNode* node = &g_game.nodes[i][j];
				if (node->used)
				{
					if (g_rotate) tsuroCardRotate(&node->card, g_rotate);
				}
			}
		}
		player->current.path = (player->current.path + g_rotate*2 + 8) % 8;
		g_rotate = 0;

		tsuroCardDump(&g_game.nodes[0][0].card);
	}
}

static tsuroNode* tsuroGetNextState(tsuroPlayer::State* nextState, const tsuroPlayer::State* state);

static void tsuroSinglePlayerUpdate(float dt)
{
	float t = g_game.player.totalLength;

	if (g_game.state == TSURO_STATE_PLAYING)
	{
		tsuroMovePlayer(dt * tsuroPlayerSpeed);

		bool stillPlaying = g_game.selectionUsed[0] || g_game.selectionUsed[1] || g_game.selectionUsed[2];

		if (stillPlaying) 
		{
			if (g_rotate)
			{
				tsuroCardRotate(&g_game.selection[g_select], g_rotate);
				tsuroSnapShot(&g_game.last);
				g_rotate = 0;
			}

			if (g_used)
			{
				// Deal a new card
				tsuroDealCard(g_select);

				// Grab undo info
				tsuroSnapShot(&g_game.last);
				g_used = 0;
			}

			if (g_use)
			{
				if (g_game.selectionUsed[g_select])
				{
					// Place the card down at the target place
					tsuroPlayer* player = &g_game.player;
					tsuroNode* node = &g_game.nodes[player->target.x][player->target.y];
					node->used = 1;
					node->card = g_game.selection[g_select];

					// Calculate the new target place
					tsuroPlayer::State nextState = player->current;
					tsuroNode* nextNode;
					while ((nextNode = tsuroGetNextState(&nextState, &nextState)) && nextNode->used) ;
					player->target = nextState;
				}
				g_use = 0;
			}
		}
		else
		{
			// Finished game by using all the tiles!
			printf("FINISHED ALL THE TILES!\n");
		}
	}
}

static void tsuroUpdate ()
{
	float dt = timerGetTickSinceLastUUpdate(&g_time);
	timerUpdate(&g_time);

	switch(g_game.mode)
	{
		case TSURO_MODE_DEMO:
			tsuroDemoUpdate(dt);
			break;

		case TSURO_MODE_SINGLE_PLAYER:
			tsuroSinglePlayerUpdate(dt);
			break;
	}

	// Reset
	g_rotate = 0;

    glutPostRedisplay();
}


static void tsuroDrawPlayer()
{
	const vec3* col = &g_game.player.colour;
	glColor3f(col->x,col->y,col->z);

	const tsuroPlayer* player = &g_game.player;
	tsuroNode* node = &g_game.nodes[player->current.x][player->current.y];
	tsuroCard* card = &node->card;

	if (node->used)
	{
		int n = (int)(player->animFraction*(float)tsuroVectorPathSize);
		if (n >= tsuroVectorPathSize) n = tsuroVectorPathSize-1;
		vec3 pos = card->vpaths[player->current.path].centre[n];
		float x = tsuroCardScreenX(player->current.x);
		float y = tsuroCardScreenY(player->current.y);
		x += pos.x;
		y += pos.y;
		drawCircle(x, y, tsuroCardSize*0.1f);
	}
	else
	{
		float x = tsuroCardScreenX(player->current.x);
		float y = tsuroCardScreenY(player->current.y);
		x += tsuroCardPoints[player->current.path][0];
		y += tsuroCardPoints[player->current.path][1];
		drawCircle(x, y, tsuroCardSize*0.1f);
	}
}
	
static void tsuroDumpPlayer()
{
	tsuroPlayer* player = &g_game.player;
	tsuroNode* node = &g_game.nodes[player->current.x][player->current.y];
	tsuroCardDump(&node->card);
}

/* Get Next Node and State from Current Position on the board */

static tsuroNode* tsuroGetNextState(tsuroPlayer::State* nextState, const tsuroPlayer::State* state)
{
	const int* other = tsuroEdgeOpposite1;

	tsuroNode* node = &g_game.nodes[state->x][state->y];
	if (node->used)
	{
		int from = state->path;
		int to = node->card.paths[state->path];

		// Transition to next node
		int destination = node->card.paths[from];
		int dx = tsuroTileTransition[destination][0];
		int dy = tsuroTileTransition[destination][1];
		int newX = state->x + dx;
		int newY = state->y + dy;
		if (newX>=0 && newY>=0 && newX<tsuroNumCells && newY<tsuroNumCells)
		{
			nextState->x = newX;
			nextState->y = newY;
			nextState->path = other[destination];
			return &g_game.nodes[nextState->x][nextState->y];
		}
		else
		{
			// Out of board, not moving anywhere
			*nextState = *state;
			return NULL;
		}
	}
	else
	{
		// Not Moving anywhere
		*nextState = *state;
		return node;
	}
}

static void tsuroMovePlayer(float frac)
{
	const int* other = tsuroEdgeOpposite1;

	tsuroPlayer* player = &g_game.player;
	tsuroNode* node = &g_game.nodes[player->current.x][player->current.y];

	if (node->used)
	{
		player->animFraction += frac;
		player->totalLength += frac;
		if (player->animFraction>=1.f)
		{
			player->animFraction -= 1.f;
			int from = player->current.path;
			int to = node->card.paths[player->current.path];

			// Colour existing path
			node->card.colour[from] = player->colour;
			node->card.colour[to] = player->colour;

			// Transition to next node
			int destination = node->card.paths[player->current.path];
			int dx = tsuroTileTransition[destination][0];
			int dy = tsuroTileTransition[destination][1];
			int newX = player->current.x + dx;
			int newY = player->current.y + dy;
			if (newX>=0 && newY>=0 && newX<tsuroNumCells && newY<tsuroNumCells)
			{
				player->current.x = newX;
				player->current.y = newY;
				player->current.path = other[destination];
			}
			else
			{
				// GAME - OVER!
				if (g_game.mode == TSURO_MODE_SINGLE_PLAYER)
				{
					player->totalLength -= player->animFraction;
					player->animFraction = 1.f;
					if (g_dragging==0)
					{
						g_game.state = TSURO_STATE_GAME_OVER;
						printf("GAME OVER!\n");
						printf("Press 'r' to reset\n");
						printf("length = %f\n", player->totalLength);
					}
				}
				if (g_game.mode == TSURO_MODE_DEMO)
				{
					tsuroReset();
				}
			}
		}
	}

	node = &g_game.nodes[player->current.x][player->current.y];
	if (!node->used)
	{
		player->totalLength -= player->animFraction;
		player->animFraction=0.f;
	}
}

static void tsuroMouseButton(int state, int x, int y)
{
	if (g_game.state == TSURO_STATE_PLAYING)
	{
		const float boardLeft = tsuroCardScreenX(0) - 0.5f*tsuroCardSize;
		const float boardRight = tsuroCardScreenX(5) + 0.5f*tsuroCardSize;
		const float boardBottom = tsuroCardScreenX(0) - 0.5f*tsuroCardSize;
		const float boardTop = tsuroCardScreenX(5) + 0.5f*tsuroCardSize;

		static int mouseDownSelection = -1;
		static int mouseDownX = -1;
		static int mouseDownY = -1;
		int selection = -1;

		// Remap y
		y = s_height - y;

		g_mouseX = x;
		g_mouseY = y;

		if (state&MOUSEDOWN)
		{
			mouseDownSelection = -1;
			mouseDownX = x;
			mouseDownY = y;
		}

		// Is the mouse being dragged
		if (!g_dragging && (mouseDownSelection>=0))
		{
			const int dragThreshold = 10;
			//if (x<boardRight)
			if (((+x-mouseDownX) > dragThreshold) || ((-x+mouseDownX) > dragThreshold) || ((+y-mouseDownY) > dragThreshold) || ((-y+mouseDownY) > dragThreshold))
			{
				// Drag card onto board and begin simulating the movement
				g_dragging = true;
			}
		}

		// Determine which card the mouse is over, if any.
		int dx = (int)(tsuroCardSize + tsuroCardSizePad);
		int x0 = (int)(tsuroCardScreenX(7) - dx * 0.5f);
		int y0 = (int)(tsuroCardScreenX(1) - dx * 0.5f);
		if ((y>=y0) && (x>=x0) && (x<(x0+dx)))
		{
			selection = (y-y0)/dx;	// 0, 1, 2
			if (selection < 3)
			{
				g_select = selection;
			}
			else
			{
				selection = -1;
			}
		}

		if (g_dragging)
		{
			if (mouseDownSelection>=0)
			{
				g_select = mouseDownSelection;
			}

			if (x>boardLeft && x<boardRight && y>boardBottom && y<boardTop)
			{
				if (!g_using)
				{
					g_use = 1;
					g_using = 1;
				}
			}
			else
			{
				g_use = 0;
				g_using = 0;
				tsuroRestore();
			}
		}

		if (state)
		{
			if (selection>=0)
			{
				if ((state==MOUSEUP) && (selection == g_select))
				{
					if (!g_dragging)
					{
						g_rotate = +1;
					}
				}
				if (state==MOUSEDOWN)
				{
					mouseDownSelection = selection;
				}
			}
			if (state==MOUSEUP)
			{
				if (g_dragging)
				{
					// Only use the card if the mouse has gone into the board area
					if (x>boardLeft && x<boardRight && y>boardBottom && y<boardTop)
					{
						g_used = true;
					}
					else
					{
						// UNDO
						tsuroRestore();
					}
				}
				g_dragging = false;
				g_using = 0;
				mouseDownSelection = -1;
			}
		}
	}
}

static void tsuroDraw()
{
	tsuroGame* game = &g_game;

    glMatrixMode (GL_MODELVIEW);
    glLoadIdentity();
	glDisable(GL_BLEND);

	float space = 0.01f*tsuroCardSize;
		
	float width=(tsuroCardSize + space)*tsuroNumCells;
	
	// Draw Grid
	if (1)
	{
		u32 N = tsuroNumCells + 1;
		float dx = tsuroCardSize + tsuroCardSizePad;
		float width = dx * (float)(N-1);
		float dy = dx;
		float x, y;
		
		glColor3f(0.2f, 0.2f, 0.2f);

		x = tsuroCardScreenX(0) - dx * 0.5f;
		y = tsuroCardScreenY(0) - dx * 0.5f;
		for (u32 i=0; i<N; i++)
		{
			vec3 from = {x, y, 0.f};
			vec3 to = {x, y+width, 0.f};
			dlDrawLine(&from, &to);
			x += dx;
		}

		x = tsuroCardScreenX(0) - dx * 0.5f;
		y = tsuroCardScreenY(0) - dx * 0.5f;
		for (u32 i=0; i<N; i++)
		{
			vec3 from = {x, y, 0.f};
			vec3 to = {x+width, y, 0.f};
			dlDrawLine(&from, &to);
			y += dy;
		}
	}

	// Draw Boundary Markers
	{
		float dx = tsuroCardSize + tsuroCardSizePad;
		float dy = dx;
		float x, y;

		// Bottom
		x = tsuroCardScreenX(0);
		y = tsuroCardScreenY(0) - dx * 0.5f;
		for (int i=0; i<tsuroNumCells; i++)
		{
			glBegin(GL_LINES);
				glVertex3f(x - (1.0/6.0)*dx, y-0.1f*dy, 0);
				glVertex3f(x - (1.0/6.0)*dx, y, 0);
				glVertex3f(x + (1.0/6.0)*dx, y-0.1f*dy, 0);
				glVertex3f(x + (1.0/6.0)*dx, y, 0);
			glEnd();
			x = x + dx;
		}

		// Top
		x = tsuroCardScreenX(0);
		y = tsuroCardScreenY(tsuroNumCells-1) + dx * 0.5f;
		for (int i=0; i<tsuroNumCells; i++)
		{
			glBegin(GL_LINES);
				glVertex3f(x - (1.0/6.0)*dx, y+0.1f*dy, 0);
				glVertex3f(x - (1.0/6.0)*dx, y, 0);
				glVertex3f(x + (1.0/6.0)*dx, y+0.1f*dy, 0);
				glVertex3f(x + (1.0/6.0)*dx, y, 0);
			glEnd();
			x = x + dx;
		}

		// Left
		x = tsuroCardScreenY(0) - dx * 0.5f;
		y = tsuroCardScreenX(0);
		for (int i=0; i<tsuroNumCells; i++)
		{
			glBegin(GL_LINES);
				glVertex3f(x-0.1f*dx, y - (1.0/6.0)*dy, 0);
				glVertex3f(x        , y - (1.0/6.0)*dy, 0);
				glVertex3f(x-0.1f*dx, y + (1.0/6.0)*dy, 0);
				glVertex3f(x        , y + (1.0/6.0)*dy, 0);
			glEnd();
			y = y + dy;
		}

		// Right
		x = tsuroCardScreenY(tsuroNumCells-1) + dx * 0.5f;
		y = tsuroCardScreenX(0);
		for (int i=0; i<tsuroNumCells; i++)
		{
			glBegin(GL_LINES);
				glVertex3f(x+0.1f*dx, y - (1.0/6.0)*dy, 0);
				glVertex3f(x        , y - (1.0/6.0)*dy, 0);
				glVertex3f(x+0.1f*dx, y + (1.0/6.0)*dy, 0);
				glVertex3f(x        , y + (1.0/6.0)*dy, 0);
			glEnd();
			y = y + dy;
		}
	}

	for (int i=0; i<tsuroNumCells; i++)
	{
		for (int j=0; j<tsuroNumCells; j++)
		{
			tsuroNode* node = &game->nodes[i][j];
			if (node->used)
			{
				tsuroCardDraw(&node->card, tsuroCardScreenX(i), tsuroCardScreenY(j));
			}
		}
	}

	tsuroDrawPlayer();

	if (g_game.mode == TSURO_MODE_SINGLE_PLAYER)
	{
		// Draw the selection board
		for (int i=0; i < 3; i++)
		{
			if (g_game.selectionUsed[i])
			{
				float x = tsuroCardScreenX(7);
				float y = tsuroCardScreenX(1+i);
				if (g_dragging && (i==g_select))
				{
					x = g_mouseX;
					y = g_mouseY;
				}
				tsuroCard* card = &g_game.selection[i];
				tsuroCardDraw(card, x, y);

				if (i == g_select)
				{
					float s = tsuroCardSize * 0.5f;
					float eps = 0.f;
					glBegin(GL_LINE_STRIP);
					glVertex3f(x+eps-s,y+eps-s,0.f);
					glVertex3f(x-eps+s,y+eps-s,0.f);
					glVertex3f(x-eps+s,y-eps+s,0.f);
					glVertex3f(x+eps-s,y-eps+s,0.f);
					glVertex3f(x+eps-s,y+eps-s,0.f);
					glEnd();
				}
			}
		}

		// Draw the simulation card
		if (g_using)
		{
			static int flash = 0;
			flash++;
			const tsuroPlayer::State* target = &g_game.last.player.target;
			float x = tsuroCardScreenX(target->x);
			float y = tsuroCardScreenX(target->y);
			float s = tsuroCardSize * 0.5f;
			float eps = 0.f;
			float alpha = float(flash%20) / 70.f;
			glColor4f(0.f,0.f,0.f,alpha);
			glEnable(GL_BLEND);
			glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
			glBegin(GL_QUADS);
			glVertex3f(x+eps-s,y+eps-s,0.f);
			glVertex3f(x-eps+s,y+eps-s,0.f);
			glVertex3f(x-eps+s,y-eps+s,0.f);
			glVertex3f(x+eps-s,y-eps+s,0.f);
			glEnd();
		}
	}
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////

static void keyboard(unsigned char key, int x, int y)
{
	if (key=='r') tsuroReset();
	if (key=='p') g_paused = !g_paused;
	if (key=='m') tsuroChangeMode();
	if (key==' ') tsuroPlayerSpeed = 4.f;
}


static void drawCircle(float x, float y, float r, bool filled)
{
    u32 N = 20;
    float dAlpha = PIx2 / (float)N;
    float a = 0.0f;
    filled ? glBegin(GL_POLYGON) : glBegin(GL_LINE_STRIP);
	for (u32 i=0; i<(N+1); i++)
	{
	    glVertex3f(x+r*cosf(a), y+r*sinf(a), 0.f);
	    a += dAlpha;
	}
    glEnd();
}

static const float drawEps = 0.003f;

static void drawPoint(float x, float y, float size = drawEps)
{
    glBegin(GL_QUADS);
	glVertex3f(x-size, y-size, 0.f);
	glVertex3f(x-size, y+size, 0.f);
	glVertex3f(x+size, y+size, 0.f);
	glVertex3f(x+size, y-size, 0.f);
    glEnd();
}

static void display()
{
    // clear the window
    glClearColor (0.5,0.5,0.5,0);
    glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // go to GL_MODELVIEW matrix mode and set the camera
    glMatrixMode (GL_MODELVIEW);
    glLoadIdentity();

    glDisable (GL_TEXTURE_2D);
    glShadeModel (GL_FLAT);
    glEnable (GL_DEPTH_TEST);
    glDepthFunc (GL_LEQUAL);
    glColor4f (1.f,1.f,1.f,1.f);

	tsuroDraw();
    
    glutSwapBuffers();
}

static void reshape(GLint width, GLint height)
{
    s_width = width;
    s_height = height;

    // setup viewport
    glViewport (0,0,s_width,s_height);
    glMatrixMode (GL_PROJECTION);
    glLoadIdentity();
	glOrtho(0.0f, width, 0.0f, height, -1.0f, 1.0f);
}

static void mouseButton(int button, int state, int x, int y)
{
	if (button==GLUT_LEFT_BUTTON)
	{
		tsuroMouseButton(state==GLUT_DOWN ? MOUSEDOWN : MOUSEUP, x, y);
	}
}

static void mouseMotion(int x, int y)
{
	tsuroMouseButton(0, x, y);
}


int main (int argc, char **argv)
{
    tsuroStart();

    // GLUT Window Initialization:
    glutInit (&argc, argv);
    glutInitWindowSize (s_width, s_height);
    glutInitDisplayMode ( GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
    glutCreateWindow ("TSURO");

    // Register callbacks:
    glutDisplayFunc (display);
    glutReshapeFunc (reshape);
    glutKeyboardFunc (keyboard);
    glutMouseFunc (mouseButton);
    glutMotionFunc (mouseMotion);
    glutPassiveMotionFunc (mouseMotion);
    glutIdleFunc (tsuroUpdate);

    //BuildPopupMenu ();
    //glutAttachMenu (GLUT_RIGHT_BUTTON);

    // Turn the flow of control over to GLUT
    glutMainLoop ();

    return 0;
}


