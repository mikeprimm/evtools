import { Component, OnInit } from '@angular/core';
import { FormControl } from '@angular/forms';
import { debounceTime, distinctUntilChanged } from 'rxjs/operators';
import { ExoplanetFetchService } from "../exoplanet-fetch.service";

@Component({
  selector: 'app-exoplanets',
  templateUrl: './exoplanets.component.html',
  styleUrls: ['./exoplanets.component.css']
})
export class ExoplanetsComponent implements OnInit  {
  public starname: FormControl;
  public debounce: number = 500;

  constructor(private api: ExoplanetFetchService) {
    this.starname = new FormControl('');
  }
  ngOnInit() {
    this.starname.valueChanges
      .pipe(debounceTime(this.debounce), distinctUntilChanged())
      .subscribe(query => {
        console.log(query);
        this.api.searchForTIC(query).subscribe(
          data => {
            console.log(data);
          },
          err => {
            console.log(err);
          }
        );
      });
  }
}
