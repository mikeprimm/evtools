import { Injectable } from '@angular/core';
import { Observable, throwError } from "rxjs";
import {
  HttpClient,
  HttpHeaders,
  HttpErrorResponse
} from "@angular/common/http";
import { catchError, map } from "rxjs/operators";

// Set the http options
const httpOptions = {
  headers: new HttpHeaders({ "Content-Type": "application/json", "Authorization": "c31z" })
};

@Injectable({
  providedIn: 'root'
})
export class ExoplanetFetchService {

  constructor(private http: HttpClient) {}

  /**
   * Function to handle error when the server return an error
   *
   * @param error
   */
  private handleError(error: HttpErrorResponse) {
    if (error.error instanceof ErrorEvent) {
      // A client-side or network error occurred. Handle it accordingly.
      console.error("An error occurred:", error.error.message);
    } else {
      // The backend returned an unsuccessful response code. The response body may contain clues as to what went wrong,
      console.error(
        `Backend returned code ${error.status}, ` + `body was: ${error.error}`
      );
    }
    // return an observable with a user-facing error message
    return throwError(error);
  }

  /**
   * Function to extract the data when the server return some
   *
   * @param res
   */
  private extractData(res: Response) {
    let body = res;
    console.log(`body=${body}`);
    return body || {};
  }

  /**
   * Function to lookup TIC from EXOFOP
   * @param searchid
   */
  public searchForTIC(searchid: string): Observable<any> {
    // Call the http GET
    return this.http.get<Response>(`https://exofop.ipac.caltech.edu/tess/gototicid.php?target=${encodeURIComponent(searchid)}&json`, httpOptions).pipe(
      map(this.extractData),
      catchError(this.handleError)
    );
  }
}
