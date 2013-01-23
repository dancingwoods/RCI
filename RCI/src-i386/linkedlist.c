#include <stdio.h>
#include <stdlib.h>
#include "linkedlist.h"

void addhead(point_t *x, point_t **head){
   if(*head != NULL){
	  x->next = *head;
	  (*head)->prev = x;
   }
   *head = x;
}


//add after point_t (not tail)
void addafter(point_t *x, point_t *a){
   x->next = a->next;
   x->prev = a;
   if(a->next != NULL)
	  a->next->prev = x;
   a->next = x;
}

void addbefore(point_t *x, point_t *b, point_t **head){
   x->next = b;
   x->prev = b->prev;
   if(b->prev != NULL)
	  b->prev->next = x;
   else
	  *head = x;
   b->prev = x;
}


//delete a point_t
void delpoint(point_t *x, point_t **head){
   if (x->prev != NULL)
	  x->prev->next = x->next;
   else
	  *head = x->next;
	  
   if (x->next != NULL)
		 x->next->prev = x->prev;
  
   free(x);
}

//removes point from list but does not free.
//always use with an add function
void removepoint(point_t *x, point_t **head){
   if (x->prev != NULL)
	  x->prev->next = x->next;
   else
	  *head = x->next;

   if (x->next != NULL)
	  x->next->prev = x->prev;
   x->next = NULL;
   x->prev = NULL;
}

void initpoint(point_t *p, double x, double y){
   p->x = x;
   p->y = y;
   p->next = NULL;
   p->prev = NULL;
}

void printlist(point_t *head){
   point_t *tmp;

   tmp = head;
   while(tmp != NULL){
	  printf("x = %g, y = %g\n", tmp->x, tmp->y);
	  tmp = tmp->next;
   }
}

/*
int main(){
   int i;
   point_t *head, *tail, *tmp;

   tmp = (point_t *)malloc(sizeof(point_t));
   initpoint(tmp, 0.0, 0.5);

   head = tail = tmp;

   for(i=1; i<10; i++){
	  tmp = (point_t *)malloc(sizeof(point_t));
	  initpoint(tmp, i, i+0.5);
	  addtail(tmp, tail);
	  tail = tmp;
   }

   printlist(head);
}
*/
